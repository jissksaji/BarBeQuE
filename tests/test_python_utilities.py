import csv
import importlib.util
import io
import os
import sys
import tempfile
import types
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"


def load_script(name, stubs=None):
    stubs = stubs or {}
    old_modules = {}
    for module_name, module in stubs.items():
        old_modules[module_name] = sys.modules.get(module_name)
        sys.modules[module_name] = module

    try:
        spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        for module_name, previous in old_modules.items():
            if previous is None:
                sys.modules.pop(module_name, None)
            else:
                sys.modules[module_name] = previous


def fake_bio_module():
    bio = types.ModuleType("Bio")
    bio.SeqIO = types.SimpleNamespace(parse=None, write=None)
    return bio


class TestParseObipcr(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parse_obipcr = load_script("parse_obipcr")

    def test_parse_fasta_header_accepts_valid_obipcr_header(self):
        parsed = self.parse_obipcr.parse_fasta_header(
            '>seq_1_sub[10..20] {"direction":"F","forward_primer":"ACGT"}'
        )

        self.assertEqual(parsed["seq_id"], "seq_1")
        self.assertEqual(parsed["start_pos"], 10)
        self.assertEqual(parsed["end_pos"], 20)
        self.assertEqual(parsed["metadata"]["direction"], "F")

    def test_parse_fasta_header_rejects_malformed_input(self):
        invalid_headers = [
            "seq_1_sub[10..20] {}",
            ">seq_1_sub[10..20]",
            ">seq_1_sub[10..20] {bad json}",
            ">seq_1[10..20] {}",
        ]

        for header in invalid_headers:
            with self.subTest(header=header):
                self.assertIsNone(self.parse_obipcr.parse_fasta_header(header))

    def test_primer_quality_helpers_cover_edge_cases(self):
        p = self.parse_obipcr

        self.assertEqual(p.calculate_gc(""), "0.0")
        self.assertEqual(p.calculate_gc("GCat"), "50.0")
        self.assertEqual(p.calculate_gc("CGAGTYTTTGAAYGCAAGTTG"), "42.86")
        self.assertEqual(p.calculate_gc("YCCCGYCTGAYCTGRGGT"), "66.67")
        self.assertEqual(p.calculate_gc("RYSWKMBDHVN"), "50.0")
        self.assertEqual(p.calculate_tm("GCGC"), "16.0")
        self.assertEqual(p.calculate_tm("GCGCGCGCGCGCGCGCGCGC"), "72.3")
        self.assertEqual(p.gc_clamp_count("AAAAAGCGCC"), "5")
        self.assertEqual(p.max_mononucleotide_run("AAGGGCTTTT"), "4")
        self.assertEqual(p.reverse_complement("acgtnx"), "NNACGT")
        self.assertEqual(p.three_prime_dimer_length("ACGT", "ACGT"), "4")
        self.assertEqual(p.max_hairpin_stem("AGCTTTTAGCT"), "4")

    def test_iupac_matching_and_mismatch_scoring(self):
        p = self.parse_obipcr

        self.assertTrue(p.is_iupac_match("Y", "C"))
        self.assertTrue(p.is_iupac_match("N", "G"))
        self.assertFalse(p.is_iupac_match("R", "C"))
        self.assertFalse(p.is_iupac_match("C", "Y"))
        self.assertFalse(p.is_iupac_match("N", "N"))

        self.assertEqual(
            p.mismatch_positions("YCC", "TYC"),
            ([2], [2]),
        )

        self.assertEqual(
            p.find_mismatches("ACGT", "ACTT", 100, is_reverse=False),
            ("3", "102", "2", "3.0"),
        )
        self.assertEqual(
            p.find_mismatches("ACGT", "ACTT", 100, is_reverse=True),
            ("3", "98", "2", "3.0"),
        )
        self.assertEqual(
            p.find_mismatches("ACNT", "ACGT", 100, is_reverse=False),
            ("", "", "", "0.0"),
        )
        self.assertEqual(
            p.find_mismatches("ACGT", "AYGT", 100, is_reverse=False),
            ("2", "101", "3", "6.0"),
        )

    def test_process_obipcr_writes_full_metrics_and_hit_spread(self):
        content = "\n".join(
            [
                '>seq1_sub[100..110] {"forward_primer":"ACGT","forward_match":"ACGT","reverse_primer":"TGCA","reverse_match":"TGCA","direction":"F","forward_error":0,"reverse_error":0}',
                "ACGTACGTACG",
                '>seq1_sub[200..214] {"forward_primer":"ACGT","forward_match":"ACCT","reverse_primer":"TGCA","reverse_match":"TGAA","direction":"R","forward_error":1,"reverse_error":1}',
                "ACGTACGTACGTACG",
                ">bad_header {}",
                "AAAA",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "obipcr.fasta"
            out = Path(tmp) / "parsed.tsv"
            inp.write_text(content + "\n")

            self.parse_obipcr.process_obipcr(inp, out)

            with out.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["Sequence_ID"], "seq1")
        self.assertEqual(rows[0]["Hits_On_Sequence"], "2")
        self.assertEqual(rows[0]["Amplicon_Length_Spread"], "4")
        self.assertEqual(rows[1]["Forward_Mismatch_Positions_Primer"], "3")
        self.assertEqual(rows[1]["Reverse_Mismatch_Positions_Genome"], "212")


class TestMask(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mask = load_script("mask")

    def test_read_fasta_joins_multiline_records(self):
        with tempfile.NamedTemporaryFile("w+", delete=False) as handle:
            handle.write(">one\nAC\nGT\n>two\nTT\n")
            path = handle.name
        try:
            self.assertEqual(
                list(self.mask.read_fasta(path)),
                [("one", "ACGT"), ("two", "TT")],
            )
        finally:
            os.unlink(path)

    def test_paired_end_masking_keeps_short_amplicons_and_masks_gap(self):
        self.assertEqual(self.mask.mask_paired_end("ACGT", 2), "ACGT")
        self.assertEqual(self.mask.mask_paired_end("AAACCCGGGTTT", 3), "AAA" + "N" * 20 + "TTT")

    def test_single_end_masking_truncates_only_long_sequences(self):
        self.assertEqual(self.mask.mask_single_end("ACGT", 10), "ACGT")
        self.assertEqual(self.mask.mask_single_end("ACGT", 2), "AC")

    def test_cli_writes_masked_fasta_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "in.fa"
            out = Path(tmp) / "out.fa"
            inp.write_text(">long\nAAACCCGGGTTT\n>short\nACGT\n")

            stdout = io.StringIO()
            with patch.object(
                sys,
                "argv",
                ["mask.py", "--input", str(inp), "--output", str(out), "--read-length", "3"],
            ), redirect_stdout(stdout):
                self.mask.main()

            self.assertEqual(out.read_text(), ">long\nAAA" + "N" * 20 + "TTT\n>short\nACGT\n")
            self.assertIn("masked   : 1", stdout.getvalue())
            self.assertIn("untouched: 1", stdout.getvalue())


class TestProcessSampleSheet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sample_sheet = load_script("process_sample_sheet")

    def test_collapse_builds_iupac_consensus(self):
        self.assertEqual(self.sample_sheet.collapse(["ACGT", "ATGT"]), "AYGT")
        self.assertEqual(self.sample_sheet.collapse(["AR", "AG"]), "AR")
        self.assertEqual(self.sample_sheet.collapse(["AZ", "AC"]), "AN")

    def test_process_sample_sheet_collapses_only_equal_length_duplicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "sheet.tsv"
            out = Path(tmp) / "collapsed.tsv"
            inp.write_text(
                "primer\tforward\treverse\tmin\tmax\n"
                "p1\tACGT\tTGCA\t10\t20\n"
                "p1\tATGT\tTGTA\t10\t20\n"
                "p2\tAAAA\tCCCC\t10\t20\n"
                "p2\tAAA\tCCC\t10\t20\n"
            )

            self.sample_sheet.process_sample_sheet(inp, out)

            self.assertEqual(
                out.read_text(),
                "primer\tforward\treverse\tmin\tmax\n"
                "p1\tAYGT\tTGYA\t10\t20\n"
                "p2\tAAAA\tCCCC\t10\t20\n"
                "p2\tAAA\tCCC\t10\t20\n",
            )


class ParsePrimersTestCase(unittest.TestCase):
    """Shared helpers for driving bin/parse_primers.py off temporary FASTA files."""

    @classmethod
    def setUpClass(cls):
        cls.parse_primers = load_script("parse_primers")

    def collect(self, files, min_len=100, max_len=500):
        """Write {name: content} into a temp dir and parse the whole directory."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        for name, content in files.items():
            (Path(tmp.name) / name).write_text(content)
        paths = self.parse_primers.resolve_inputs(Path(tmp.name))
        return self.parse_primers.collect_rows(paths, min_len, max_len)


class TestParsePrimersNaming(ParsePrimersTestCase):
    def test_single_pair_file_is_named_after_the_file(self):
        rows, warnings, errors = self.collect({"ITS2.fasta": ">ITS2_fwd\nACGT\n>ITS2_rev\nTGCA\n"})

        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])
        self.assertEqual(rows, [{"primer": "ITS2", "fwd": "ACGT", "rev": "TGCA", "min": 100, "max": 500}])

    def test_two_prefixes_stay_separate_even_when_lengths_match(self):
        rows, _warnings, errors = self.collect(
            {"markers.fasta": ">MA_FWD\nAAAA\n>MA_REV\nTTTT\n>POL_FWD\nCCCC\n>POL_REV\nGGGG\n"}
        )

        self.assertEqual(errors, [])
        self.assertEqual(
            [(r["primer"], r["fwd"], r["rev"]) for r in rows],
            [("markers_MA", "AAAA", "TTTT"), ("markers_POL", "CCCC", "GGGG")],
        )

    def test_prefix_and_number_combine_when_one_prefix_splits(self):
        rows, _warnings, errors = self.collect(
            {
                "markers.fasta": ">MA_FWD_1\nAAAA\n>MA_FWD_2\nAAAAAA\n>MA_REV\nTTTT\n"
                ">POL_FWD\nCCCC\n>POL_REV\nGGGG\n"
            }
        )

        self.assertEqual(errors, [])
        self.assertEqual(
            [r["primer"] for r in rows],
            ["markers_MA_1", "markers_MA_2", "markers_POL"],
        )


class TestParsePrimersCollapsing(ParsePrimersTestCase):
    def test_same_length_variants_collapse_into_one_degenerate_primer(self):
        rows, warnings, errors = self.collect({"MA.fasta": ">MA_fwd_1\nACGT\n>MA_fwd_2\nATGT\n>MA_rev\nTGCA\n"})

        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])
        self.assertEqual([(r["primer"], r["fwd"], r["rev"]) for r in rows], [("MA", "AYGT", "TGCA")])

    def test_differing_lengths_split_into_numbered_sets_and_warn(self):
        rows, warnings, errors = self.collect({"ITS2.fasta": ">ITS2_fwd_1\nAAAA\n>ITS2_fwd_2\nAAAAAA\n>ITS2_rev\nTTTT\n"})

        self.assertEqual(errors, [])
        self.assertEqual(
            [(r["primer"], r["fwd"], r["rev"]) for r in rows],
            [("ITS2_1", "AAAA", "TTTT"), ("ITS2_2", "AAAAAA", "TTTT")],
        )
        self.assertEqual(len(warnings), 1)
        self.assertIn("ITS2.fasta", warnings[0])
        self.assertIn("[4, 6]", warnings[0])

    def test_both_directions_splitting_gives_every_combination(self):
        rows, warnings, _errors = self.collect(
            {"X.fasta": ">X_fwd_1\nAAAA\n>X_fwd_2\nAAAAAA\n>X_rev_1\nTTTT\n>X_rev_2\nTTTTTT\n"}
        )

        self.assertEqual(
            [(r["primer"], r["fwd"], r["rev"]) for r in rows],
            [
                ("X_1", "AAAA", "TTTT"),
                ("X_2", "AAAA", "TTTTTT"),
                ("X_3", "AAAAAA", "TTTT"),
                ("X_4", "AAAAAA", "TTTTTT"),
            ],
        )
        self.assertEqual(len(warnings), 1)

    def test_two_record_file_without_direction_tokens_is_fwd_then_rev(self):
        rows, _warnings, errors = self.collect({"pair.fasta": ">first\nAAAA\n>second\nTTTT\n"})

        self.assertEqual(errors, [])
        self.assertEqual([(r["primer"], r["fwd"], r["rev"]) for r in rows], [("pair", "AAAA", "TTTT")])


class TestParsePrimersRejects(ParsePrimersTestCase):
    def assert_single_error(self, files, *fragments):
        rows, _warnings, errors = self.collect(files)

        self.assertEqual(rows, [])
        self.assertEqual(len(errors), 1)
        for fragment in fragments:
            self.assertIn(fragment, errors[0])

    def test_rejects_empty_sequence(self):
        self.assert_single_error({"bad.fasta": ">bad_fwd\n\n>bad_rev\nTTTT\n"}, "bad.fasta", "bad_fwd", "empty")

    def test_rejects_non_nucleotide_characters(self):
        self.assert_single_error({"bad.fasta": ">bad_fwd\nACGX\n>bad_rev\nTTTT\n"}, "bad.fasta", "bad_fwd", "X")

    def test_rejects_file_with_no_fasta_records(self):
        self.assert_single_error({"bad.fasta": "ACGTACGT\n"}, "bad.fasta", "no FASTA records")

    def test_rejects_empty_file(self):
        self.assert_single_error({"bad.fasta": ""}, "bad.fasta", "no FASTA records")

    def test_rejects_prefix_missing_a_direction(self):
        self.assert_single_error(
            {"bad.fasta": ">bad_fwd_1\nAAAA\n>bad_fwd_2\nCCCC\n>bad_fwd_3\nGGGG\n"},
            "bad.fasta",
            "3 fwd",
            "0 rev",
        )

    def test_rejects_untagged_record_mixed_with_tagged_records(self):
        self.assert_single_error(
            {"bad.fasta": ">p_fwd\nAAAA\n>p_rev\nTTTT\n>extra\nGGGG\n"}, "bad.fasta", "extra"
        )

    def test_collects_errors_from_every_file_before_giving_up(self):
        rows, _warnings, errors = self.collect(
            {
                "good.fasta": ">good_fwd\nAAAA\n>good_rev\nTTTT\n",
                "bad_one.fasta": ">x_fwd\nACGX\n>x_rev\nTTTT\n",
                "bad_two.fasta": ">y_fwd\nAAAA\n",
            }
        )

        self.assertEqual(rows, [])
        self.assertEqual(len(errors), 2)
        self.assertTrue(any("bad_one.fasta" in e for e in errors))
        self.assertTrue(any("bad_two.fasta" in e for e in errors))


class TestParsePrimersInputs(ParsePrimersTestCase):
    def test_directory_picks_up_fasta_extensions_and_ignores_everything_else(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ["a.fasta", "b.fa", "c.fna", "reads.fastq.gz", "notes.txt"]:
                (Path(tmp) / name).write_text(">p_fwd\nAAAA\n>p_rev\nTTTT\n")

            paths = self.parse_primers.resolve_inputs(Path(tmp))

            self.assertEqual([p.name for p in paths], ["a.fasta", "b.fa", "c.fna"])

    def test_single_fasta_file_is_accepted_directly(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "ITS2.fasta"
            fasta.write_text(">ITS2_fwd\nACGT\n>ITS2_rev\nTGCA\n")

            paths = self.parse_primers.resolve_inputs(fasta)
            rows, _warnings, errors = self.parse_primers.collect_rows(paths, 100, 500)

            self.assertEqual(paths, [fasta])
            self.assertEqual(errors, [])
            self.assertEqual([r["primer"] for r in rows], ["ITS2"])

    def test_directory_without_any_fasta_files_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "reads.fastq.gz").write_text("")

            with self.assertRaises(SystemExit):
                self.parse_primers.resolve_inputs(Path(tmp))

    def test_is_fasta_distinguishes_fasta_from_a_samplesheet(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "primers.fasta"
            sheet = Path(tmp) / "sheet.tsv"
            fasta.write_text("\n\n>p_fwd\nAAAA\n")
            sheet.write_text("primer\tfwd\trev\tmin\tmax\n")

            self.assertTrue(self.parse_primers.is_fasta(fasta))
            self.assertFalse(self.parse_primers.is_fasta(sheet))

    def test_amplicon_bounds_are_written_onto_every_row(self):
        rows, _warnings, errors = self.collect(
            {"a.fasta": ">a_fwd\nAAAA\n>a_rev\nTTTT\n", "b.fasta": ">b_fwd\nCCCC\n>b_rev\nGGGG\n"},
            min_len=120,
            max_len=460,
        )

        self.assertEqual(errors, [])
        self.assertEqual([(r["min"], r["max"]) for r in rows], [(120, 460), (120, 460)])


class TestParsePrimersMain(ParsePrimersTestCase):
    def test_main_writes_a_samplesheet_and_a_warnings_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "ITS2.fasta").write_text(">ITS2_fwd_1\nAAAA\n>ITS2_fwd_2\nAAAAAA\n>ITS2_rev\nTTTT\n")
            out = Path(tmp) / "primers.tsv"
            warnings = Path(tmp) / "primer_warnings.txt"

            with patch.object(
                sys,
                "argv",
                [
                    "parse_primers.py",
                    "--input", tmp,
                    "--min", "100",
                    "--max", "500",
                    "--out", str(out),
                    "--warnings", str(warnings),
                ],
            ), redirect_stderr(io.StringIO()):
                self.parse_primers.main()

            self.assertEqual(
                out.read_text(),
                "primer\tfwd\trev\tmin\tmax\n"
                "ITS2_1\tAAAA\tTTTT\t100\t500\n"
                "ITS2_2\tAAAAAA\tTTTT\t100\t500\n",
            )
            self.assertIn("ITS2.fasta", warnings.read_text())

    def test_main_exits_and_writes_nothing_when_a_file_is_invalid(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "bad.fasta").write_text(">bad_fwd\nACGX\n>bad_rev\nTTTT\n")
            out = Path(tmp) / "primers.tsv"

            with patch.object(
                sys,
                "argv",
                ["parse_primers.py", "--input", tmp, "--min", "100", "--max", "500", "--out", str(out)],
            ):
                with self.assertRaises(SystemExit):
                    self.parse_primers.main()

            self.assertFalse(out.exists())


class TestTaxidAndAccessionHelpers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.taxid_filter = load_script("taxid_db_filter", {"Bio": fake_bio_module()})

    def test_taxid_db_filter_loads_requested_taxids_and_matching_accessions(self):
        with tempfile.TemporaryDirectory() as tmp:
            taxids = Path(tmp) / "taxids.txt"
            mapping = Path(tmp) / "accession_taxid.tsv"
            taxids.write_text("\n111\n222\n333\n\n")
            mapping.write_text("A1.1\t111\nA2\t999\nA3.12\t222\nAY846379.1.1791\t333\na malformed line\n")

            keep_taxids = self.taxid_filter.load_taxids(taxids)
            accessions = self.taxid_filter.load_matching_accessions(mapping, keep_taxids)

        self.assertEqual(keep_taxids, {"111", "222", "333"})
        self.assertEqual(accessions, {"A1", "A3", "AY846379"})


class TestAccessionBlocklist(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.accession_filter = load_script("filter_accession_blocklist")

    def test_loads_comments_and_normalizes_accession_versions(self):
        with tempfile.NamedTemporaryFile("w+", delete=False, suffix=".txt") as handle:
            handle.write(
                "# accessions to ignore\n"
                "mk123456.1  # versioned and lower-case\n"
                "\n"
                "AY846379.1.1791\n"
                "MK123456\n"
            )
            path = handle.name
        try:
            accessions = self.accession_filter.load_accession_blocklist(path)
        finally:
            os.unlink(path)

        self.assertEqual(accessions, {"MK123456", "AY846379"})

    def test_filters_parsed_tsv_and_fasta_and_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            blocklist = tmp / "accession_blocklist.txt"
            fasta = tmp / "amplicons.fasta"
            parsed = tmp / "parsed.tsv"
            fasta_output = tmp / "filtered.fasta"
            tsv_output = tmp / "filtered.tsv"
            summary = tmp / "summary.tsv"

            blocklist.write_text("# remove these\nA1\na3.1\nNOT_PRESENT\n")
            fasta.write_text(
                ">A1.1 blocked by base accession\nAAAA\n"
                ">A2.1 kept\nCC\nCC\n"
                ">A3.5 blocked case-insensitively\nGGGG\n"
            )
            parsed.write_text(
                "Sequence_ID\tAmplicon_Length\n"
                "A1.1\t4\n"
                "A2.1\t4\n"
                "A3.5\t4\n"
                "A3.5\t5\n"
            )

            blocked = self.accession_filter.load_accession_blocklist(blocklist)
            fasta_total, fasta_removed, fasta_matched = (
                self.accession_filter.filter_fasta(fasta, fasta_output, blocked)
            )
            tsv_total, tsv_removed, tsv_matched = self.accession_filter.filter_tsv(
                parsed, tsv_output, blocked
            )
            self.accession_filter.write_summary(
                summary,
                blocked,
                (fasta_total, fasta_removed),
                (tsv_total, tsv_removed),
                fasta_matched | tsv_matched,
            )

            self.assertEqual(fasta_output.read_text(), ">A2.1 kept\nCC\nCC\n")
            self.assertEqual(
                tsv_output.read_text(),
                "Sequence_ID\tAmplicon_Length\nA2.1\t4\n",
            )
            summary_text = summary.read_text()
            self.assertIn("listed_accessions\t3", summary_text)
            self.assertIn("matched_accessions\t2", summary_text)
            self.assertIn("unmatched_accessions\t1", summary_text)
            self.assertIn("fasta_removed_records\t2", summary_text)
            self.assertIn("parsed_removed_rows\t3", summary_text)

    def test_main_rejects_an_empty_accession_blocklist(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            blocklist = tmp / "accession_blocklist.txt"
            fasta = tmp / "amplicons.fasta"
            parsed = tmp / "parsed.tsv"
            blocklist.write_text("# no accessions\n\n")
            fasta.write_text(">A1.1\nAAAA\n")
            parsed.write_text("Sequence_ID\tAmplicon_Length\nA1.1\t4\n")

            with patch.object(
                sys,
                "argv",
                [
                    "filter_accession_blocklist.py",
                    "--fasta", str(fasta),
                    "--tsv", str(parsed),
                    "--accession-blocklist", str(blocklist),
                    "--fasta-output", str(tmp / "out.fasta"),
                    "--tsv-output", str(tmp / "out.tsv"),
                    "--summary", str(tmp / "summary.tsv"),
                ],
            ):
                with self.assertRaises(SystemExit):
                    self.accession_filter.main()

    def test_rejects_a_parsed_table_without_sequence_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            parsed = tmp / "parsed.tsv"
            parsed.write_text("accession\tvalue\nA1\t4\n")
            with self.assertRaisesRegex(ValueError, "Sequence_ID"):
                self.accession_filter.filter_tsv(
                    parsed,
                    tmp / "out.tsv",
                    {"A1"},
                )


class TestDbDistribution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fake_taxid_tools = types.ModuleType("taxidTools")
        cls.db_distribution = load_script("db_distribution", {"taxidTools": fake_taxid_tools})

    def test_read_taxid_counts_skips_headers_and_bad_rows(self):
        with tempfile.NamedTemporaryFile("w+", delete=False) as handle:
            handle.write("taxid\tcount\n123\t4\nbad\t5\n456\tseven\n789\t1\n")
            path = handle.name
        try:
            self.assertEqual(self.db_distribution.read_taxid_counts(path), {123: 4, 789: 1})
        finally:
            os.unlink(path)

    def test_resolve_lineages_walks_to_root_and_handles_missing_taxids(self):
        root = types.SimpleNamespace(taxid="1", rank="no rank", name="root", parent=None)
        genus = types.SimpleNamespace(taxid="10", rank="genus", name="Genus", parent=root)
        species = types.SimpleNamespace(taxid="11", rank="species", name="Genus species", parent=genus)
        tax = {"1": root, "10": genus, "11": species}

        lineages, ranks, names = self.db_distribution.resolve_lineages(tax, [11, 999])

        self.assertEqual(lineages[11], [1, 10, 11])
        self.assertEqual(lineages[999], [])
        self.assertEqual(ranks[11], "species")
        self.assertEqual(names[10], "Genus")

    def test_write_tsv_outputs_ranked_lineage_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "distribution.tsv"
            self.db_distribution.write_tsv(
                {11: 3},
                {11: [1, 10, 11]},
                {1: "no rank", 10: "genus", 11: "species"},
                {1: "root", 10: "Genus", 11: "Genus species"},
                out,
            )

            self.assertEqual(
                out.read_text(),
                "taxid\tcount\tresolved_rank\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"
                "11\t3\tspecies\t\t\t\t\t\tGenus\tGenus species\n",
            )


class TestCompletenessTable(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fake_ete3 = types.ModuleType("ete3")
        fake_ete3.NCBITaxa = object
        cls.completeness = load_script("completeness_table", {"ete3": fake_ete3})

    def test_load_species_list_from_inline_string_or_file(self):
        self.assertEqual(
            self.completeness.load_species_list("Camellia sinensis + Camellia japonica"),
            ["Camellia sinensis", "Camellia japonica"],
        )
        with tempfile.NamedTemporaryFile("w+", delete=False) as handle:
            handle.write("A species\n\nB species\n")
            path = handle.name
        try:
            self.assertEqual(self.completeness.load_species_list(path), ["A species", "B species"])
        finally:
            os.unlink(path)

    def test_build_row_uses_lineage_and_db_counts(self):
        class FakeNcbi:
            lineage = {
                111: [1, 10, 100, 111],
                112: [1, 10, 100, 112],
                211: [1, 20, 200, 211],
            }
            ranks = {
                1: "no rank",
                10: "family",
                20: "family",
                100: "genus",
                200: "genus",
                111: "species",
                112: "species",
                211: "species",
            }
            names = {
                10: "Theaceae",
                100: "Camellia",
                111: "Camellia sinensis",
                112: "Camellia japonica",
                211: "Other species",
            }

            def get_name_translator(self, names):
                return {"Camellia sinensis": [111]}

            def get_lineage(self, taxid):
                return self.lineage[taxid]

            def get_rank(self, taxids):
                return {taxid: self.ranks[taxid] for taxid in taxids}

            def get_descendant_taxa(self, taxid, collapse_subspecies=False):
                if taxid == 100:
                    return [111, 112]
                if taxid == 10:
                    return [111, 112, 211]
                return []

            def get_taxid_translator(self, taxids):
                return {taxid: self.names[taxid] for taxid in taxids}

        row = self.completeness.build_row("Camellia sinensis", {111: 7, 112: 2}, FakeNcbi())

        self.assertEqual(
            row,
            [
                "Camellia sinensis",
                7,
                2,
                2,
                100.0,
                2,
                3,
                66.67,
                "Camellia japonica; Camellia sinensis",
                "Camellia japonica; Camellia sinensis",
            ],
        )


if __name__ == "__main__":
    unittest.main()
