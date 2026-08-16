#!/usr/bin/env python3
"""Export the FULL database records for the Camellia sinensis ITS2_collapsed hits.

The amplicon export (export_camellia_its2_amplicons.py) contains only the region
obipcr cut out. This script pulls the complete source sequences for the same
accessions straight out of the reference FASTA the pipeline ran against, and
joins them to the taxonomy / amplicon coordinates.

Usage:
    python export_camellia_its2_full_sequences.py [--db euphyllophyta.fasta]
                                                  [--amplicons camellia_..._amplicons.tsv]
                                                  [--fasta extracted.fasta]
                                                  [--out out.tsv]

--fasta is a cache: if the file exists it is reused, otherwise the records are
extracted from --db with `seqkit grep` (a full scan of a ~60 GB file).
"""
import argparse
import io
import subprocess
from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).parent.parent / "results"
DEFAULT_DB = Path("/home/saj/jiss/ncbi_extraction/euphyllophyta.fasta")

parser = argparse.ArgumentParser()
parser.add_argument("--db", type=Path, default=DEFAULT_DB)
parser.add_argument("--amplicons", type=Path, default=RESULTS / "camellia_sinensis_ITS2_collapsed_amplicons.tsv")
parser.add_argument("--fasta", type=Path, default=RESULTS / "camellia_sinensis_ITS2_collapsed_full_sequences.fasta")
parser.add_argument("--out", type=Path, default=RESULTS / "camellia_sinensis_ITS2_collapsed_full_sequences.tsv")
args = parser.parse_args()

amplicons = pd.read_csv(args.amplicons, sep="\t", dtype=str)

if not args.fasta.exists():
    ids_file = args.fasta.with_suffix(".ids.txt")
    ids_file.write_text("\n".join(amplicons["Sequence_ID"]) + "\n")
    print(f"extracting {len(amplicons)} records from {args.db} (full scan, slow) ...")
    with args.fasta.open("w") as handle:
        subprocess.run(["seqkit", "grep", "-f", str(ids_file), str(args.db)], stdout=handle, check=True)

# seqkit fx2tab gives us id / description / sequence / length without hand-parsing FASTA.
fx2tab = subprocess.run(
    ["seqkit", "fx2tab", "--name", "--only-id", "--length", str(args.fasta)],
    capture_output=True, text=True, check=True,
)
lengths = pd.read_csv(io.StringIO(fx2tab.stdout), sep="\t", names=["Sequence_ID", "Full_Sequence_Length"], dtype=str)

sequences = subprocess.run(
    ["seqkit", "fx2tab", str(args.fasta)],
    capture_output=True, text=True, check=True,
)
records = pd.read_csv(
    io.StringIO(sequences.stdout),
    sep="\t", names=["header", "Full_Sequence", "qual"], dtype=str,
)
records["Sequence_ID"] = records["header"].str.split(n=1).str[0]
records["Definition"] = records["header"].str.split(n=1).str[1]
records = records.merge(lengths, on="Sequence_ID")

merged = amplicons.merge(records, on="Sequence_ID", how="left", validate="one_to_one")

columns = [
    "accession", "Sequence_ID", "accession_taxid", "accession_name",
    "cluster_id", "assigned_name", "assigned_taxid", "assigned_rank",
    "Definition", "Full_Sequence_Length",
    "Forward_Binding_Start", "Reverse_Binding_End", "Amplicon_Length", "Direction",
    "Amplicon_Sequence", "Full_Sequence",
]
merged = merged[columns].sort_values(["accession_name", "accession"])
merged.to_csv(args.out, sep="\t", index=False)

missing = merged["Full_Sequence"].isna().sum()
print(f"wrote {args.out} — {len(merged)} records, {missing} without a database sequence")
print(merged["accession_name"].value_counts().to_string())
print(merged["Full_Sequence_Length"].astype(float).describe()[["min", "50%", "max"]].to_string())
