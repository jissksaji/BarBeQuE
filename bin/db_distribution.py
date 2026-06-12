#!/usr/bin/env python3
"""
db_distribution.py — BarBeQuE DB composition table
Resolves taxids from the reference DB to their full lineage.

Usage:
    db_distribution.py --input <taxids_counts.tsv> --output <out.tsv>
"""

import argparse
from ete3 import NCBITaxa

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


# I/O

def read_taxid_counts(path):
    """Read taxid\tcount TSV, return {taxid: count}."""

    db_taxids = {}

    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2 or parts[0].lower() == "taxid":
                continue
            try:
                db_taxids[int(parts[0])] = int(parts[1])
            except ValueError:
                continue

    return db_taxids


def write_tsv(db_taxids, all_lineages, all_ranks, all_names, output_path):
    """Write one row per taxid with full lineage and sequence count."""

    with open(output_path, "w") as out:

        out.write("taxid\tcount\tresolved_rank\t" + "\t".join(RANKS) + "\n")

        for taxid, count in sorted(db_taxids.items()):

            lin = all_lineages.get(taxid, [])

            ranked_lin = [
                (tid, all_ranks.get(tid, ""), all_names.get(tid, "?"))
                for tid in lin
                if all_ranks.get(tid) in RANKS
            ]

            rank_map = {rank: name for _, rank, name in ranked_lin}
            my_rank  = all_ranks.get(taxid, "no rank")

            cols = (
                [str(taxid), str(count), my_rank]
                + [rank_map.get(r, "") for r in RANKS]
            )

            out.write("\t".join(cols) + "\n")


# Taxonomy

def resolve_lineages(ncbi, taxid_list):
    """
    Batch-resolve lineages for all taxids.
    Returns all_lineages, all_ranks, all_names.
    """

    all_lineages = {}
    for taxid in taxid_list:
        try:
            all_lineages[taxid] = ncbi.get_lineage(taxid)
        except Exception:
            all_lineages[taxid] = []

    all_tids = set()
    for lin in all_lineages.values():
        all_tids.update(lin)

    all_ranks = ncbi.get_rank(list(all_tids))
    all_names = ncbi.get_taxid_translator(list(all_tids))

    return all_lineages, all_ranks, all_names


# Main

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True,  help="taxid counts TSV from BUILD_DB_TAXIDS")
    parser.add_argument("--output", required=True,  help="output TSV path")
    parser.add_argument("--dbfile", default=None,   help="ete3 sqlite path (optional)")
    args = parser.parse_args()

    ncbi = NCBITaxa(dbfile=args.dbfile) if args.dbfile else NCBITaxa()

    db_taxids = read_taxid_counts(args.input)
    all_lineages, all_ranks, all_names = resolve_lineages(ncbi, list(db_taxids.keys()))
    write_tsv(db_taxids, all_lineages, all_ranks, all_names, args.output)


if __name__ == "__main__":
    main()