import argparse
import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "modules" / "helper" / "cluster_consensus" / "cluster_consensus.py"


def load_cluster_consensus():
    taxid_tools = types.ModuleType("taxidTools")
    previous = sys.modules.get("taxidTools")
    sys.modules["taxidTools"] = taxid_tools
    try:
        spec = importlib.util.spec_from_file_location("cluster_consensus", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if previous is None:
            sys.modules.pop("taxidTools", None)
        else:
            sys.modules["taxidTools"] = previous


class FakeNode:
    def __init__(self, name="Species A", taxid="101", rank="species"):
        self.name = name
        self.taxid = taxid
        self.rank = rank


class RecordingTaxonomy:
    def __init__(self, result=None):
        self.names = {
            "101": "Species A",
            "202": "Species B",
            "missing": None,
        }
        self.result = result if result is not None else FakeNode()
        self.calls = []

    def getName(self, taxid):
        return self.names.get(taxid)

    def consensus(self, taxids, **kwargs):
        self.calls.append((taxids, kwargs))
        return self.result


class TestClusterConsensus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_cluster_consensus()

    def test_consensus_fraction_accepts_supported_range(self):
        self.assertEqual(self.module.consensus_fraction("0.8"), 0.8)
        self.assertEqual(self.module.consensus_fraction("1.0"), 1.0)

    def test_consensus_fraction_rejects_ambiguous_or_invalid_values(self):
        for value in ("0.5", "0.4", "1.01"):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    self.module.consensus_fraction(value)

    def test_fraction_and_accession_votes_are_passed_to_taxidtools(self):
        tax = RecordingTaxonomy()
        taxids = ["101"] * 8 + ["202"] * 2 + ["missing"]

        assignment = self.module.cluster_consensus(tax, taxids, 0.8)

        self.assertEqual(
            assignment,
            ["Species A", "101", "species", "Species A;Species B"],
        )
        self.assertEqual(len(tax.calls), 1)
        voted_taxids, options = tax.calls[0]
        self.assertEqual(voted_taxids, ["101"] * 8 + ["202"] * 2)
        self.assertEqual(options["min_consensus"], 0.8)
        self.assertTrue(options["ignore_missing"])

    def test_default_fraction_preserves_strict_consensus(self):
        args = self.module.parse_args(
            ["--input", "input.tsv", "--taxdump", "taxdump", "--output", "out.tsv"]
        )
        self.assertEqual(args.min_consensus, 1.0)

    def test_cluster_without_valid_taxids_is_unclassified(self):
        tax = RecordingTaxonomy()

        assignment = self.module.cluster_consensus(tax, ["missing"], 0.8)

        self.assertEqual(
            assignment,
            ["Unclassified", "Unknown", "no rank", ""],
        )
        self.assertEqual(tax.calls, [])

    def test_missing_consensus_node_is_unclassified(self):
        tax = RecordingTaxonomy(result=False)
        tax.result = None

        assignment = self.module.cluster_consensus(tax, ["101", "202"], 0.8)

        self.assertEqual(
            assignment,
            ["Unclassified", "Unknown", "no rank", "Species A;Species B"],
        )


if __name__ == "__main__":
    unittest.main()
