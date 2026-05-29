#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#ete.py
#This script does tree traversal.
# traverses from family to species.
# any node in the tree not present in the blast db is set as not present
# pass means present in the blast db and found
# fail means present in the blast db and not found
# no_data means not present in the blast db
#the tree use ete to get the current no of entries in each rank 

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

    # We pass a taxon name, and must translate it to an NCBI id
    name2taxid = ncbi.get_name_translator([taxname])

    # the first entry, of which the first value is taken
    try:
        taxid = list(name2taxid.items())[0][1][0]
    except IndexError:
        message = f"{taxname} is not a valid taxa name in the NCBI nomenclature. Common errors include misspelling and use of a taxid or common name instead of the latin name."
        raise ValueError(message)

    # the consensus report
    r = open(report, "r")
    lines = r.readlines()
    bucket = {}

    for line in lines:
        line = line.strip()
        if not line:
            continue
        elements = line.split("\t")
        if len(elements) >= 5:  
            sci = elements[-4]
            tax = elements[-3]
            bucket[tax] = sci

    r.close()

    # the list taxids included in the BlastDB
    with open(refs, "r") as taxids:
        blast_tax = taxids.readlines()
        blast_tax = [tax.rstrip("\n") for tax in blast_tax]

    ok = "#7ee076"
    fail = "#dfc2b1"
    missing = "#eeeeee"
    tree = ncbi.get_descendant_taxa(taxid, collapse_subspecies=False, return_tree=True)
    species_names = [
        n.sci_name for n in tree.traverse()
        if n.rank == "species" and re.match(r"^[A-Z][a-z]*\s[a-z]*$", n.sci_name)
    ]
    species_count = len(species_names)
    print(f"Total species in tree: {species_count}", flush=True)
    with open(output + ".species_count.txt", "w") as sc:
        sc.write(f"root_taxon\t{taxname}\n")
        sc.write(f"total_species_in_tree\t{species_count}\n")
        for name in sorted(species_names):
            sc.write(f"{name}\n")
            
    # Calculate Family Coverage before modifying tree nodes
    family_data = {}
    for n in tree.traverse():
        if n.rank == "family":
            family_tid = str(n.name)
            descendant_tids = [str(desc.name) for desc in n.traverse()]
            if any(tid in bucket for tid in descendant_tids):
                family_data[n.sci_name] = f"OK\t{family_tid}\t{ok}\tfamily"
            elif any(tid in blast_tax for tid in descendant_tids):
                family_data[n.sci_name] = f"FAIL\t{family_tid}\t{fail}\tfamily"
            else:
                family_data[n.sci_name] = f"NO_DATA\t{family_tid}\t{missing}\tfamily"

    # Calculate Genus Coverage
    genus_data = {}
    for n in tree.traverse():
        if n.rank == "genus":
            genus_tid = str(n.name)
            descendant_tids = [str(desc.name) for desc in n.traverse()]
            if any(tid in bucket for tid in descendant_tids):
                genus_data[n.sci_name] = f"OK\t{genus_tid}\t{ok}\tgenus"
            elif any(tid in blast_tax for tid in descendant_tids):
                genus_data[n.sci_name] = f"FAIL\t{genus_tid}\t{fail}\tgenus"
            else:
                genus_data[n.sci_name] = f"NO_DATA\t{genus_tid}\t{missing}\tgenus"

    data = {}

    nodes = []
    for n in tree.traverse():
        if n.rank == "species":
            tid = str(n.name)

            # NCBI taxonomy is full of non-species level terminal leafs
            # we skip all leafs not matching the 'Genus species' pattern
            if re.match(r"^[A-Z][a-z]*\s[a-z]*$", n.sci_name):
                # This taxon is in the blast db and was found
                if tid in bucket:
                    data[n.sci_name] = f"OK\t{tid}\t{ok}\tspecies"
                # This taxon is in the blast db and was not found
                elif tid in blast_tax:
                    data[n.sci_name] = f"FAIL\t{tid}\t{fail}\tspecies"
                # this taxon was not in the blast db
                else:
                    data[n.sci_name] = f"NO_DATA\t{tid}\t{missing}\tspecies"
                
                n.name = n.sci_name
                nodes.append(n)

    tree.prune(nodes)

    print(tree.write(outfile=output + ".nwk", format=1))

    f = open(output + ".tsv", "w")
    f.write("Taxon\tStatus\tTaxid\tColor\tRank\n")
    for taxon in family_data:
        status = family_data[taxon]
        f.write(f"{taxon}\t{status}\n")
    for taxon in genus_data:
        status = genus_data[taxon]
        f.write(f"{taxon}\t{status}\n")
    for taxon in data:
        status = data[taxon]
        f.write(f"{taxon}\t{status}\n")

    f.close()


if __name__ == "__main__":
    main(args.taxon, args.reference, args.report, args.output)
