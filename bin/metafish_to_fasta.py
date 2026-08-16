#!/usr/bin/env python3

"""Convert the downloaded MetaFish CSV file to FASTA."""

import csv
import re
import sys


def clean_sequence(sequence):
    sequence = sequence or ""
    return re.sub(r"[^ACGTURYSWKMBDHVNacgturyswkmbdhvn]", "", sequence).upper().replace("U", "T")


input_csv, output_fasta = sys.argv[1:3]

with open(input_csv, newline="") as source, open(output_fasta, "w") as output:
    for row in csv.DictReader(source):
        accession = (row.get("gbAccession") or row.get("dbid") or "").strip()
        sequence = clean_sequence(row.get("nucleotides"))

        if not accession or not sequence:
            continue

        species = (row.get("sciNameValid") or row.get("sciNameOrig") or "").strip()
        header = f"{accession} {species}".strip()

        output.write(f">{header}\n")
        for start in range(0, len(sequence), 80):
            output.write(sequence[start : start + 80] + "\n")
