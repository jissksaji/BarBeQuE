#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Legacy strict-LCA cluster consensus implementation.

This file preserves the implementation used before fractional consensus was
introduced. It is not called by the CLUSTER_CONSENSUS Nextflow process.
"""

import os
import argparse
import taxidTools


parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--input", help="Path to cluster_taxonomy.tsv from JOIN_ACCESSION_TAXONOMY")
parser.add_argument("--taxdump", help="Path to the taxdump folder")
parser.add_argument("--output", help="Path to output table")
args = parser.parse_args()


def load_taxonomy(taxdump):
    return taxidTools.read_taxdump(
        os.path.join(taxdump, "nodes.dmp"),
        os.path.join(taxdump, "rankedlineage.dmp"),
        os.path.join(taxdump, "merged.dmp")
    )


def parse_cluster_taxonomy(input_file):
    """Parse cluster_taxonomy.tsv (cluster_id, accession, taxid, phylum, class, order, family, genus, species)"""
    clusters_taxid = {}    # cluster_id -> list of taxids
    rows = []              # keep original rows for output

    with open(input_file, 'r') as fi:
        for line in fi.readlines():
            elements = line.rstrip('\n').split('\t')
            cluster_id = elements[0]
            accession = elements[1]
            taxid = elements[2]

            clusters_taxid.setdefault(cluster_id, [])
            clusters_taxid[cluster_id].append(taxid)

            rows.append(elements)

    return clusters_taxid, rows


def main(input_file, taxdump, output):
    tax = load_taxonomy(taxdump)

    # Parse input
    clusters_taxid, rows = parse_cluster_taxonomy(input_file)

    # Get consensus, rank, and disambiguation per cluster
    clusters_cons = {}
    for cluster_id, taxids in clusters_taxid.items():
        # Get unique taxids and their names for disambiguation
        unique_taxids = list(dict.fromkeys(taxids))
        disambiguation = ";".join([name for name in [tax.getName(t) for t in unique_taxids] if name is not None])
        valid_taxids = [t for t in taxids if tax.getName(t) is not None]

        if valid_taxids:
            cons = tax.lca(valid_taxids, ignore_missing=True)
            clusters_cons[cluster_id] = [
                cons.name,
                str(cons.taxid),
                cons.rank,
                disambiguation,
            ]
        else:
            clusters_cons[cluster_id] = [
                "Unclassified",
                "Unknown",
                "no rank",
                disambiguation,
            ]

    # Dump: original row + accession_name + LCA columns appended
    with open(output, 'w') as fo:
        for elements in rows:
            cluster_id = elements[0]
            accession_taxid = elements[2]
            node = tax.get(accession_taxid)
            accession_name = node.name if node else "Unknown"

            line = "\t".join(elements + [accession_name] + clusters_cons[cluster_id])
            fo.write(line + "\n")


if __name__ == "__main__":
    main(
        args.input,
        args.taxdump,
        args.output
    )
