#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from ete3 import NCBITaxa
import re
import sys

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--taxon", help="Taxonomy to use as root")
parser.add_argument("--reference", help="Taxonomy lists from the DB")
parser.add_argument("--report", help="The consensus report")
parser.add_argument("--output")
args = parser.parse_args()

ncbi = NCBITaxa()


def main(taxname, refs, report, output):

    # Accept either a numeric taxid or a Latin name
    # ete2 acccepts anything as long as its valid , dosen't matter the rank
    if str(taxname).isdigit():
        taxid = int(taxname)
    else:
        name2taxid = ncbi.get_name_translator([taxname])
        try:
            taxid = list(name2taxid.items())[0][1][0]
        except IndexError:
            raise ValueError(
                f"{taxname} is not a valid taxa name in the NCBI nomenclature. "
                "Common errors include misspelling and use of a taxid or common name instead of the latin name."
            )

    # the consensus report: taxid -> accession
    bucket = {}
    with open(report, "r") as r:
        for line in r:
            elements = line.split("\t")
            sci = elements[1]   # Accession
            tax = elements[2]   # Taxid
            bucket[tax] = sci

    # the list taxids included in the BlastDB
    with open(refs, "r") as taxids:
        blast_tax = taxids.readlines()
        blast_tax = [tax.rstrip('\n') for tax in blast_tax]

    ok = "#7ee076"
    fail = "#dfc2b1"
    missing = "#eeeeee"
    tree = ncbi.get_descendant_taxa(taxid, collapse_subspecies=False, return_tree=True,intermediate_nodes=True)
    data = {}

    nodes = []
    for n in tree.traverse(strategy="postorder"):
        tid = str(n.name)
        
        # Determine presence in consensus report (bucket) or blast DB
        n.is_in_bucket = (tid in bucket) or any(getattr(c, "is_in_bucket", False) for c in n.children)
        n.is_in_db = (tid in blast_tax) or any(getattr(c, "is_in_db", False) for c in n.children)
        
        # Only process standard taxonomic ranks
        if n.rank in ["species", "subspecies", "varietas", "forma", "genus", "family", "order", "class", "phylum", "kingdom"]:
            if n.is_in_bucket:
                data[n.sci_name] = f"OK\t{tid}\t{ok}\t{n.rank}"
            elif n.is_in_db:
                data[n.sci_name] = f"FAIL\t{tid}\t{fail}\t{n.rank}"
            else:
                data[n.sci_name] = f"NO_DATA\t{tid}\t{missing}\t{n.rank}"
            nodes.append(n)
                
        # Rename node to scientific name for the .nwk tree output
        n.name = n.sci_name

    tree.prune(nodes)
    
    # Add missing leaf nodes for taxa in bucket but not in tree 
    for tax, sci in bucket.items():
        if sci not in data:
            rank_dict = ncbi.get_rank([int(tax)])
            rank = rank_dict.get(int(tax), "no rank")
            data[sci] = f"OK\t{tax}\t{ok}\t{rank}"

    tree.write(outfile=output + ".nwk", format=1, format_root_node=True)

    f = open(output + ".tsv", "w")
    f.write("Taxon\tStatus\tTaxid\tColor\tRank\n")
    for taxon in data:
        status = data[taxon]
        f.write(f"{taxon}\t{status}\n")

    f.close()


if __name__ == '__main__':
    main(args.taxon, args.reference, args.report, args.output)
