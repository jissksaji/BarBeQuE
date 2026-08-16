import csv
import importlib.util
import io
import os
import sys
import tempfile
import types
import unittest
from contextlib import redirect_stdout
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


class TestCollapsePrimers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collapse_primers = load_script("collapse_primers")

    def test_classify_detects_direction_tokens_case_insensitively(self):
        c = self.collapse_primers

        self.assertEqual(c.classify("primer_FWD_1"), "fwd")
        self.assertEqual(c.classify("primer reverse 1"), "rev")
        self.assertIsNone(c.classify("primer1"))

    def test_main_collapses_tagged_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "primers.fa"
            out = Path(tmp) / "out.fa"
            fasta.write_text(">p_forward_a\nACGT\n>p_forward_b\nATGT\n>p_reverse_a\nTGCA\n>p_reverse_b\nTGTA\n")

            with patch.object(
                sys,
                "argv",
                ["collapse_primers.py", "--fasta", str(fasta), "--prefix", "P", "--out", str(out)],
            ):
                self.collapse_primers.main()

            self.assertEqual(out.read_text(), ">P_fwd\nAYGT\n>P_rev\nTGYA\n")

    def test_main_accepts_plain_two_record_fasta_without_collapsing(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "primers.fa"
            out = Path(tmp) / "out.fa"
            fasta.write_text(">first\nAAAA\n>second\nTTTT\n")

            with patch.object(
                sys,
                "argv",
                ["collapse_primers.py", "--fasta", str(fasta), "--prefix", "P", "--out", str(out)],
            ):
                self.collapse_primers.main()

            self.assertEqual(out.read_text(), ">P_fwd\nAAAA\n>P_rev\nTTTT\n")

    def test_main_rejects_ambiguous_or_unequal_primer_sets(self):
        with tempfile.TemporaryDirectory() as tmp:
            ambiguous = Path(tmp) / "ambiguous.fa"
            unequal = Path(tmp) / "unequal.fa"
            out = Path(tmp) / "out.fa"
            ambiguous.write_text(">one\nAAAA\n>two\nTTTT\n>three\nCCCC\n")
            unequal.write_text(">one_fwd\nAAAA\n>two_fwd\nAAA\n>one_rev\nTTTT\n")

            for fasta in [ambiguous, unequal]:
                with self.subTest(fasta=fasta.name), patch.object(
                    sys,
                    "argv",
                    ["collapse_primers.py", "--fasta", str(fasta), "--prefix", "P", "--out", str(out)],
                ):
                    with self.assertRaises(SystemExit):
                        self.collapse_primers.main()


class TestTaxidAndAccessionHelpers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.taxid_filter = load_script("taxid_db_filter", {"Bio": fake_bio_module()})
        cls.accession_filter = load_script("filter_accessions", {"Bio": fake_bio_module()})

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

    def test_filter_accessions_ignores_comments_and_strips_versions(self):
        with tempfile.NamedTemporaryFile("w+", delete=False) as handle:
            handle.write("# comment\nMK123456.1\n\nAB2\n")
            path = handle.name
        try:
            self.assertEqual(self.accession_filter.load_exclusions(path), {"MK123456", "AB2"})
        finally:
            os.unlink(path)


class TestBlocklistFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.blocklist_filter = load_script("filter_blocklist")

    def test_filters_fasta_from_foodme_taxids_and_writes_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            blocklist = tmp / "blocklist.txt"
            mapping = tmp / "accession_taxid.tsv"
            fasta = tmp / "database.fasta"
            output = tmp / "filtered.fasta"
            summary = tmp / "summary.tsv"

            blocklist.write_text("# FooDMe2 blocklist\n111 # unwanted\n333\n")
            mapping.write_text(
                "accession\taccession.version\ttaxid\tgi\n"
                "A1\tA1.1\t111\t0\n"
                "A2\tA2.1\t222\t0\n"
                "NOT_IN_FASTA\tNOT_IN_FASTA.1\t111\t0\n"
                "A3.5\t333\n"
            )
            fasta.write_text(
                ">A1.1 blocked four-column mapping\nAAAA\n"
                ">A2.1 kept\nCC\nCC\n"
                ">A3.5 blocked two-column mapping\nGGGG\n"
            )

            blocked_taxids = self.blocklist_filter.load_blocked_taxids(blocklist)
            database_accessions = self.blocklist_filter.load_fasta_accessions(fasta)
            blocked_accessions = self.blocklist_filter.load_blocked_accessions(
                mapping, blocked_taxids, database_accessions
            )
            total, removed = self.blocklist_filter.filter_fasta(
                fasta, output, blocked_accessions
            )
            self.blocklist_filter.write_summary(
                summary, total, removed, blocked_taxids, blocked_accessions
            )

            self.assertEqual(blocked_taxids, {"111", "333"})
            self.assertEqual(database_accessions, {"A1", "A2", "A3"})
            self.assertEqual(blocked_accessions, {"A1", "A3"})
            self.assertEqual(output.read_text(), ">A2.1 kept\nCC\nCC\n")
            self.assertIn("kept_records\t1", summary.read_text())
            self.assertIn("removed_records\t2", summary.read_text())


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
