import csv
import importlib.util
import io
import json
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

    def test_gc_calculation_handles_iupac_bases(self):
        p = self.parse_obipcr

        self.assertEqual(p.calculate_gc(""), "0.0")
        self.assertEqual(p.calculate_gc("GCat"), "50.0")
        self.assertEqual(p.calculate_gc("CGAGTYTTTGAAYGCAAGTTG"), "42.86")
        self.assertEqual(p.calculate_gc("YCCCGYCTGAYCTGRGGT"), "66.67")
        self.assertEqual(p.calculate_gc("RYSWKMBDHVN"), "50.0")

    def test_iupac_matching_and_mismatch_positions(self):
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
        self.assertEqual(rows[0]["Amplicon_GC"], "54.55")
        self.assertEqual(rows[1]["Forward_Mismatch_Positions_Primer"], "3")
        self.assertEqual(rows[1]["Reverse_Mismatch_From_3Prime"], "2")


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
        paths = self.parse_primers.find_input_files(Path(tmp.name))
        return self.parse_primers.collect_rows(paths, min_len, max_len)


class TestParsePrimersNaming(ParsePrimersTestCase):
    def test_single_pair_file_is_named_after_the_file(self):
        rows, warnings, errors = self.collect({"ITS2.fasta": ">ITS2_fwd\nACGT\n>ITS2_rev\nTGCA\n"})

        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])
        self.assertEqual(rows, [{"primer": "ITS2", "fwd": "ACGT", "rev": "TGCA", "min": 100, "max": 500}])

    def test_one_file_is_always_exactly_one_primer_set(self):
        rows, _warnings, errors = self.collect(
            {"markers.fasta": ">MA_FWD\nAAAA\n>MA_REV\nTTTT\n>POL_FWD\nCCCC\n>POL_REV\nGGGG\n"}
        )

        self.assertEqual(errors, [])
        self.assertEqual(
            [(r["primer"], r["fwd"], r["rev"]) for r in rows],
            [("markers", "MMMM", "KKKK")],
        )

    def test_variant_tag_before_the_direction_is_the_same_primer(self):
        # FooDMe2's 16S_ASU184.fasta names its alternative reverse MA_ALT_REV,
        # i.e. the variant tag sits between the name and the direction.
        rows, _warnings, errors = self.collect(
            {
                "16S_ASU184.fasta": ">MA_FWD\nGACGAGAAGACCCTATGGAGC\n"
                ">MA_REV\nTCCGAGGTCACCCCAACC\n"
                ">POL_FWD\nGACGAGAAGACCCTGTGGAAC\n"
                ">POL_REV\nTCCAAGGTCGCCCCAACC\n"
                ">MA_ALT_REV\nTCCGAGATCACCCCAATC\n"
            }
        )

        self.assertEqual(errors, [])
        self.assertEqual(
            [(r["primer"], r["fwd"], r["rev"]) for r in rows],
            [("16S_ASU184", "GACGAGAAGACCCTRTGGARC", "TCCRAGRTCRCCCCAAYC")],
        )


class TestParsePrimersCollapsing(ParsePrimersTestCase):
    def test_same_length_variants_collapse_into_one_degenerate_primer(self):
        rows, warnings, errors = self.collect({"MA.fasta": ">MA_fwd_1\nACGT\n>MA_fwd_2\nATGT\n>MA_rev\nTGCA\n"})

        self.assertEqual(errors, [])
        self.assertEqual([(r["primer"], r["fwd"], r["rev"]) for r in rows], [("MA", "AYGT", "TGCA")])
        # The samplesheet now holds a sequence that is in no input file, so say so.
        self.assertEqual(len(warnings), 1)
        self.assertIn("MA.fasta", warnings[0])

    def test_a_single_pair_collapses_without_a_warning(self):
        rows, warnings, errors = self.collect({"MA.fasta": ">MA_fwd\nACGT\n>MA_rev\nTGCA\n"})

        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])
        self.assertEqual([(r["fwd"], r["rev"]) for r in rows], [("ACGT", "TGCA")])

    def test_differing_fwd_lengths_are_rejected(self):
        rows, _warnings, errors = self.collect(
            {"ITS2.fasta": ">ITS2_fwd_1\nAAAA\n>ITS2_fwd_2\nAAAAAA\n>ITS2_rev\nTTTT\n"}
        )

        self.assertEqual(rows, [])
        self.assertEqual(len(errors), 1)
        self.assertIn("ITS2.fasta", errors[0])
        self.assertIn("fwd primers have different lengths", errors[0])
        self.assertIn("4, 6", errors[0])

    def test_differing_rev_lengths_are_rejected(self):
        rows, _warnings, errors = self.collect(
            {"X.fasta": ">X_fwd\nAAAA\n>X_rev_1\nTTTT\n>X_rev_2\nTTTTTT\n"}
        )

        self.assertEqual(rows, [])
        self.assertEqual(len(errors), 1)
        self.assertIn("rev primers have different lengths", errors[0])

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

    def test_rejects_file_missing_a_direction(self):
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

            paths = self.parse_primers.find_input_files(Path(tmp))

            self.assertEqual([p.name for p in paths], ["a.fasta", "b.fa", "c.fna"])

    def test_single_fasta_file_is_accepted_directly(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "ITS2.fasta"
            fasta.write_text(">ITS2_fwd\nACGT\n>ITS2_rev\nTGCA\n")

            paths = self.parse_primers.find_input_files(fasta)
            rows, _warnings, errors = self.parse_primers.collect_rows(paths, 100, 500)

            self.assertEqual(paths, [fasta])
            self.assertEqual(errors, [])
            self.assertEqual([r["primer"] for r in rows], ["ITS2"])

    def test_several_paths_can_be_passed_at_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "folder"
            folder.mkdir()
            (folder / "c.fasta").write_text(">c_fwd\nAAAA\n>c_rev\nTTTT\n")
            first = Path(tmp) / "a.fasta"
            second = Path(tmp) / "b.fasta"
            first.write_text(">a_fwd\nAAAA\n>a_rev\nTTTT\n")
            second.write_text(">b_fwd\nCCCC\n>b_rev\nGGGG\n")

            paths = self.parse_primers.find_input_files([first, second, folder])
            rows, _warnings, errors = self.parse_primers.collect_rows(paths, 100, 500)

            self.assertEqual(errors, [])
            self.assertEqual([r["primer"] for r in rows], ["a", "b", "c"])

    def test_two_files_with_the_same_name_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            for folder in ("one", "two"):
                (Path(tmp) / folder).mkdir()
                (Path(tmp) / folder / "ITS2.fasta").write_text(">p_fwd\nAAAA\n>p_rev\nTTTT\n")

            paths = self.parse_primers.find_input_files(
                [Path(tmp) / "one", Path(tmp) / "two"]
            )
            rows, _warnings, errors = self.parse_primers.collect_rows(paths, 100, 500)

            self.assertEqual(rows, [])
            self.assertEqual(len(errors), 1)
            self.assertIn("ITS2", errors[0])

    def test_directory_without_any_fasta_files_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "reads.fastq.gz").write_text("")

            with self.assertRaises(SystemExit):
                self.parse_primers.find_input_files(Path(tmp))

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
            (Path(tmp) / "ITS2.fasta").write_text(">ITS2_fwd_1\nAAAA\n>ITS2_fwd_2\nACAA\n>ITS2_rev\nTTTT\n")
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
                "ITS2\tAMAA\tTTTT\t100\t500\n",
            )
            self.assertIn("ITS2.fasta", warnings.read_text())

    def test_main_accepts_several_input_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "a.fasta").write_text(">a_fwd\nAAAA\n>a_rev\nTTTT\n")
            (Path(tmp) / "b.fasta").write_text(">b_fwd\nCCCC\n>b_rev\nGGGG\n")
            out = Path(tmp) / "primers.tsv"

            with patch.object(
                sys,
                "argv",
                [
                    "parse_primers.py",
                    "--input", str(Path(tmp) / "a.fasta"), str(Path(tmp) / "b.fasta"),
                    "--min", "100",
                    "--max", "500",
                    "--out", str(out),
                ],
            ), redirect_stdout(io.StringIO()):
                self.parse_primers.main()

            self.assertEqual(
                out.read_text(),
                "primer\tfwd\trev\tmin\tmax\n"
                "a\tAAAA\tTTTT\t100\t500\n"
                "b\tCCCC\tGGGG\t100\t500\n",
            )

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


# Two real variants of the same length. They combine to ACSTW, which also
# accepts ACGTT and ACCTA - sequences no real primer has.
VARIANT_PRIMERS = """>v_fwd
ACGTA
>v_rev
ACCTT
"""

# Only one 5 nt primer, so a 5 nt binding site is decided by that primer alone.
SINGLE_LENGTH_PRIMERS = """>only_fwd
ACGTA
>only_rev
GGGGGGGG
"""


class FilterObipcrTestCase(unittest.TestCase):
    """Shared helpers for driving bin/filter_obipcr.py off temporary files."""

    @classmethod
    def setUpClass(cls):
        cls.filter_obipcr = load_script("filter_obipcr")

    def amplicon_header(self, name, fwd_primer, fwd_match, rev_primer, rev_match):
        """One raw obipcr header, in the shape obipcr actually writes."""
        annotation = {
            "definition": "test record",
            "direction": "forward",
            "forward_error": 0,
            "forward_match": fwd_match,
            "forward_primer": fwd_primer,
            "reverse_error": 0,
            "reverse_match": rev_match,
            "reverse_primer": rev_primer,
        }
        return f">{name}_sub[1..10] {json.dumps(annotation)}"

    def run_filter(self, primers, headers, mismatches=0, primer_name="MARK"):
        """Run the CLI over one primer FASTA and return (kept headers, stats)."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        directory = Path(tmp.name)

        primer_dir = directory / "primers"
        primer_dir.mkdir()
        (primer_dir / f"{primer_name}.fasta").write_text(primers)

        record_text = ""
        for header in headers:
            record_text += f"{header}\nGATTACAGAT\n"

        amplicons = directory / "raw.fasta"
        amplicons.write_text(record_text)

        out = directory / "filtered.fasta"
        stats = directory / "stats.tsv"

        argv = [
            "filter_obipcr.py",
            "--amplicons", str(amplicons),
            "--primer-input", str(primer_dir),
            "--primer", primer_name,
            "--db", "testdb",
            "--mismatches", str(mismatches),
            "--out", str(out),
            "--stats", str(stats),
        ]
        with patch.object(sys, "argv", argv), redirect_stderr(io.StringIO()):
            self.filter_obipcr.main()

        kept = []
        for line in out.read_text().splitlines():
            if line.startswith(">"):
                kept.append(line)

        return kept, stats.read_text()


class TestFilterObipcrMatching(FilterObipcrTestCase):
    def test_ambiguous_target_base_counts_as_a_mismatch(self):
        find = self.filter_obipcr.find_mismatch_positions
        # An unknown reference base is no evidence that a primer would bind, and
        # obipcr counts it as a mismatch - verified against its own error counts.
        self.assertEqual(find("ACGTA", "ACGTN"), [4])
        # Even when the primer carries the very same ambiguity code.
        self.assertEqual(find("ACGTY", "ACGTY"), [4])
        # A degenerate primer still matches a concrete base that it covers.
        self.assertEqual(find("ACGTY", "ACGTC"), [])

    def test_clamped_positions_come_from_the_obipcr_pattern(self):
        self.assertEqual(
            self.filter_obipcr.get_clamped_positions("ACGT#C#"),
            ("ACGTC", {3, 4}),
        )

    def test_matching_ignores_the_direction_label(self):
        # Both variants are 5 nt, so the rev-labelled primer is available to
        # explain a forward binding site.
        self.assertTrue(
            self.filter_obipcr.is_site_explained(
                "ACCTT", ["ACGTA", "ACCTT"], 0, set(), "MARK"
            )
        )

    def test_a_site_is_only_judged_against_primers_of_its_own_length(self):
        # The ITS2 shape: 21 nt forward variant, 18 nt reverse variant. An 18 nt
        # site is explained by the 18 nt primer; the 21 nt one is simply skipped.
        its2 = ["CGAGTYTTTGAAYGCAAGTTG", "YCCCGYCTGAYCTGRGGT"]
        self.assertTrue(
            self.filter_obipcr.is_site_explained(
                "CCCCGCCTGACCTGAGGT", its2, 0, set(), "ITS2"
            )
        )

    def test_a_site_with_no_primer_of_its_length_is_an_error(self):
        # No primer of this length means the wrong primer FASTA was paired with
        # this obipcr output, which must not look like a working filter.
        with self.assertRaises(SystemExit) as raised:
            self.filter_obipcr.is_site_explained(
                "ACGT", ["ACGTA", "ACGTAC"], 0, set(), "MARK"
            )
        self.assertIn("does not belong", str(raised.exception))


class TestFilterObipcrRecords(FilterObipcrTestCase):
    def test_a_site_only_the_combined_primer_explains_is_discarded(self):
        chimera = self.amplicon_header("chimera", "ACSTW", "acgtt", "ACSTW", "acgta")
        real = self.amplicon_header("real", "ACSTW", "acgta", "ACSTW", "acctt")

        kept, _ = self.run_filter(VARIANT_PRIMERS, [chimera, real])

        # Only the record whose both sites are real variants survives, and its
        # header is written back unchanged for PARSE_OBIPCR.
        self.assertEqual(kept, [real])

    def test_a_mismatch_on_a_clamped_base_discards_within_the_budget(self):
        clamped = self.amplicon_header(
            "clamped", "ACGT#A#", "acgtt", "GGGGGGGG", "gggggggg"
        )
        free = self.amplicon_header(
            "free", "ACGTA", "acgtt", "GGGGGGGG", "gggggggg"
        )

        kept, _ = self.run_filter(SINGLE_LENGTH_PRIMERS, [clamped], mismatches=1)
        self.assertEqual(kept, [])

        # The same single mismatch passes once the clamp is gone.
        kept, _ = self.run_filter(SINGLE_LENGTH_PRIMERS, [free], mismatches=1)
        self.assertEqual(kept, [free])

    def test_stats_file_records_the_discards(self):
        chimera = self.amplicon_header("chimera", "ACSTW", "acgtt", "ACSTW", "acgta")
        real = self.amplicon_header("real", "ACSTW", "acgta", "ACSTW", "acctt")

        _, stats = self.run_filter(VARIANT_PRIMERS, [chimera, real])

        rows = stats.strip().splitlines()
        self.assertEqual(
            rows[0].split("\t"),
            [
                "primer",
                "db",
                "total",
                "kept",
                "discarded",
                "forward_unexplained",
                "reverse_unexplained",
            ],
        )
        self.assertEqual(
            rows[1].split("\t"),
            ["MARK", "testdb", "2", "1", "1", "1", "0"],
        )

    def test_a_primer_set_with_no_matching_fasta_is_reported(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        primer_dir = Path(tmp.name)
        (primer_dir / "MARK.fasta").write_text(VARIANT_PRIMERS)

        with self.assertRaises(SystemExit) as raised:
            self.filter_obipcr.find_primer_fasta(primer_dir, "OTHER")

        self.assertIn("OTHER", str(raised.exception))
        self.assertIn("MARK", str(raised.exception))

    def test_the_primer_fasta_is_found_by_its_file_name(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        primer_dir = Path(tmp.name)
        (primer_dir / "ITS2.fasta").write_text(VARIANT_PRIMERS)
        (primer_dir / "trnL.fasta").write_text(VARIANT_PRIMERS)

        # PARSE_PRIMERS names each set after its file, so meta.primer selects it.
        found = self.filter_obipcr.find_primer_fasta(primer_dir, "trnL")
        self.assertEqual(found.name, "trnL.fasta")


class TestTaxidAndAccessionHelpers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.taxid_filter = load_script("taxid_db_filter", {"Bio": fake_bio_module()})

    def test_taxid_db_filter_loads_requested_taxids_and_matching_accessions(self):
        with tempfile.TemporaryDirectory() as tmp:
            taxids = Path(tmp) / "taxids.txt"
            mapping = Path(tmp) / "accession_taxid.tsv"
            taxids.write_text("\n111\n222\n333\n\n")
            mapping.write_text(
                "A1.1\t111\n"
                "A2\t999\n"
                "A3.12\t222\n"
                "AY846379.1.1791\t333\n"
                "NCBI1\tNCBI1.4\t111\t0\n"
                "NCBI2\tNCBI2.1\t999\t0\n"
                "a malformed line\n"
            )

            keep_taxids = self.taxid_filter.load_taxids(taxids)
            accessions = self.taxid_filter.load_matching_accessions(mapping, keep_taxids)

        self.assertEqual(keep_taxids, {"111", "222", "333"})
        self.assertEqual(accessions, {"A1", "A3", "AY846379", "NCBI1"})


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


if __name__ == "__main__":
    unittest.main()
