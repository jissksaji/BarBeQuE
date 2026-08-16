#!/usr/bin/env python3
"""Classify every teeliste hit by the kind of GenBank record it was cut from.

For each of the 21 primers in the teeliste matrix, this finds the hits belonging
to the 83 teeliste species (species + infraspecific taxa, matched by name prefix)
and looks up the source record's definition line in the obipcr raw FASTA headers,
which carry the original GenBank definition as JSON.

The point is to see how many hits come from genome assemblies rather than from
targeted barcode submissions, and for which species the assembly-derived records
were the only thing that produced a hit at all.

Usage: python analyze_teeliste_hit_sources.py [results_dir]
"""
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

RESULTS = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent.parent / "results")
DB = "euphyllophyta"

CONSENSUS_COLUMNS = [
    "cluster_id", "accession", "accession_taxid", "accession_name",
    "assigned_name", "assigned_taxid", "assigned_rank", "disambiguation",
]

# A definition counts as assembly-derived only on genome-level wording. Mentioning
# "chloroplast" is not enough: targeted plastid barcodes ("rbcL ... partial cds;
# chloroplast") name the organelle too, and those are ordinary submissions.
GENOME_PATTERN = (r"(complete|draft|partial|whole)[^,;]{0,30}genome|genome assembly|"
                  r"\bchromosome\b|\bscaffold\b|\bcontig\b|plastome")
PREDICTED_PATTERN = r"^predicted:|, mrna$|, ncrna$|transcript variant|uncharacterized loc"
BARCODE_PATTERN = (r"internal transcribed spacer|ribosomal rna|rrna|its\d|genomic dna containing|"
                   r"maturase|matk|rbcl|trn|psb|rpo|atp|spacer|gene[s]?[ ,]|cds")

primers = [p.strip() for p in
           pd.read_csv(RESULTS / "teeliste_83_by_21_primers_matrix.tsv", sep="\t", nrows=0).columns[5:]]

teeliste = pd.read_csv(RESULTS / "teeliste.tsv", sep="\t", dtype=str)
teeliste.columns = ["german_name", "latin_name", "taxid", "rank", "basis"]
species_names = sorted(teeliste["latin_name"].dropna().unique(), key=len, reverse=True)
name_pattern = re.compile("^(" + "|".join(re.escape(n) for n in species_names) + ")\\b")


def species_of(accession_name):
    """Map an accession's taxon name onto its teeliste species (varieties included)."""
    match = name_pattern.match(accession_name or "")
    return match.group(1) if match else None


# 1. Collect the teeliste hits per primer from the cluster consensus files.
hits = []
for primer in primers:
    path = RESULTS / "consensus" / f"{primer}_{DB}.cluster_consensus.tsv"
    if not path.exists():
        print(f"  skipping {primer}: no consensus file")
        continue
    df = pd.read_csv(path, sep="\t", names=CONSENSUS_COLUMNS, dtype=str).drop_duplicates()
    df["species"] = df["accession_name"].map(species_of)
    df = df[df["species"].notna()].copy()
    df["primer"] = primer
    hits.append(df[["primer", "species", "accession", "accession_taxid", "accession_name", "cluster_id"]])

hits = pd.concat(hits, ignore_index=True)
wanted = set(hits["accession"])
print(f"{len(hits)} teeliste hits across {hits['primer'].nunique()} primers, {len(wanted)} distinct accessions")

# 2. Pull the GenBank definition for those accessions out of the raw obipcr headers.
# Scanning ~2.5 GB of headers is slow, so reuse definitions from a previous run.
definitions = {}
cache = RESULTS / "teeliste_hit_sources.tsv"
if cache.exists():
    cached = pd.read_csv(cache, sep="\t", dtype=str).dropna(subset=["definition"])
    definitions = dict(zip(cached["accession"], cached["definition"]))
    print(f"reused {len(definitions)} cached definitions from {cache.name}")

header_re = re.compile(r"^>(\S+?)(?:_sub\[[^\]]*\])?\s+(\{.*\})\s*$")
for primer in primers if wanted - set(definitions) else []:
    path = RESULTS / "raw" / "obipcr" / f"{primer}_{DB}_raw.fasta"
    if not path.exists():
        continue
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            match = header_re.match(line)
            if not match:
                continue
            accession = match.group(1).split(".")[0]
            if accession in wanted and accession not in definitions:
                try:
                    definitions[accession] = json.loads(match.group(2)).get("definition", "")
                except json.JSONDecodeError:
                    pass

hits["definition"] = hits["accession"].map(definitions)
print(f"definitions resolved for {hits['definition'].notna().sum()} / {len(hits)} hits")


def categorise(definition):
    if not isinstance(definition, str) or not definition:
        return "unknown"
    text = definition.lower()
    if re.search(GENOME_PATTERN, text):
        if re.search(r"chloroplast|plastid", text):
            return "plastid_genome"
        if re.search(r"mitochondri", text):
            return "mitochondrial_genome"
        return "nuclear_assembly"
    if re.search(PREDICTED_PATTERN, text):
        return "predicted_from_assembly"
    if re.search(BARCODE_PATTERN, text):
        return "targeted_barcode"
    return "other"


hits["source_category"] = hits["definition"].map(categorise)
hits["from_assembly"] = hits["source_category"].isin(
    ["plastid_genome", "mitochondrial_genome", "nuclear_assembly", "predicted_from_assembly"])

hits.to_csv(RESULTS / "teeliste_hit_sources.tsv", sep="\t", index=False)

print("\n=== hits by source category ===")
print(hits["source_category"].value_counts().to_string())
print(f"\nassembly-derived hits: {hits['from_assembly'].sum()} / {len(hits)} "
      f"({100 * hits['from_assembly'].mean():.1f}%)")

# 3. Where was it critical? Species whose only hits come from assembly records.
per_species = hits.groupby("species").agg(
    hits=("accession", "size"),
    accessions=("accession", "nunique"),
    primers=("primer", "nunique"),
    assembly_hits=("from_assembly", "sum"),
).reset_index()
per_species["assembly_only"] = per_species["assembly_hits"] == per_species["hits"]
per_species.to_csv(RESULTS / "teeliste_hit_sources_by_species.tsv", sep="\t", index=False)

detected = set(per_species["species"])
missing = sorted(set(teeliste["latin_name"]) - detected)
print(f"\nteeliste species with at least one hit: {len(detected)} / {len(teeliste)}")
if missing:
    print("no hits at all: " + ", ".join(missing))

critical = per_species[per_species["assembly_only"]]
print(f"\n=== species detectable ONLY via assembly-derived records: {len(critical)} ===")
print(critical[["species", "hits", "accessions", "primers"]].to_string(index=False))

# Per primer, how much of the signal is assembly-derived.
print("\n=== assembly share per primer ===")
per_primer = hits.groupby("primer").agg(hits=("accession", "size"), assembly=("from_assembly", "sum"))
per_primer["pct"] = (100 * per_primer["assembly"] / per_primer["hits"]).round(1)
print(per_primer.sort_values("pct", ascending=False).to_string())
