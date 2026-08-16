#!/usr/bin/env python3
"""Export the ITS2_collapsed in-silico amplicons for Camellia sinensis.

Covers the species itself (taxid 4442) plus its infraspecific varieties, joining
the cluster consensus taxonomy to the parsed obipcr amplicons.

Usage: python export_camellia_its2_amplicons.py [results_dir] [output.tsv]
"""
import sys
from pathlib import Path

import pandas as pd

STEM = "ITS2_collapsed_euphyllophyta"
SPECIES_PREFIX = "Camellia sinensis"

CONSENSUS_COLUMNS = [
    "cluster_id",
    "accession",
    "accession_taxid",
    "accession_name",
    "assigned_name",
    "assigned_taxid",
    "assigned_rank",
    "disambiguation",
]

results_dir = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent.parent / "results")
out_path = Path(sys.argv[2] if len(sys.argv) > 2 else results_dir / "camellia_sinensis_ITS2_collapsed_amplicons.tsv")

consensus = pd.read_csv(
    results_dir / "consensus" / f"{STEM}.cluster_consensus.tsv",
    sep="\t",
    names=CONSENSUS_COLUMNS,
    dtype=str,
).drop_duplicates()

camellia = consensus[consensus["accession_name"].str.startswith(SPECIES_PREFIX, na=False)]

amplicons = pd.read_csv(results_dir / "parsed_obipcr" / f"{STEM}.tsv", sep="\t", dtype=str)
# obipcr keeps the GenBank version suffix (KY928292.1); the taxonomy join does not.
amplicons["accession"] = amplicons["Sequence_ID"].str.split(".").str[0]

merged = camellia.merge(amplicons, on="accession", how="left")

columns = [
    "accession",
    "Sequence_ID",
    "accession_taxid",
    "accession_name",
    "cluster_id",
    "assigned_name",
    "assigned_taxid",
    "assigned_rank",
    "Amplicon_Length",
    "Forward_Binding_Start",
    "Reverse_Binding_End",
    "Direction",
    "Forward_Primer",
    "Reverse_Primer",
    "Forward_Errors",
    "Reverse_Errors",
    "Amplicon_GC",
    "Hits_On_Sequence",
    "Amplicon_Sequence",
]
merged = merged[columns].sort_values(["accession_name", "accession"])
merged.to_csv(out_path, sep="\t", index=False)

print(f"wrote {out_path} — {len(merged)} amplicon rows, {merged['accession'].nunique()} accessions")
print(merged["accession_name"].value_counts().to_string())
print(f"missing amplicon sequence: {merged['Amplicon_Sequence'].isna().sum()}")
