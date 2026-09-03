#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys

def get_args():
    parser = argparse.ArgumentParser(description="Augment taxonomic coverage table with representation metrics.")
    parser.add_argument("--coverage", required=True, help="Path to coverage.tsv (output of ete.py)")
    parser.add_argument("--consensus", required=True, help="Path to consensus.tsv (output of CLUSTER_CONSENSUS)")
    parser.add_argument("--output", required=True, help="Path to output representation.tsv")
    return parser.parse_args()

def main():
    args = get_args()

    # 2. Single pass over consensus.tsv
    aggregates = {}

    with open(args.consensus, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 7:
                continue # skip rows with <7 columns
            
            # col2: taxid (index 2)
            # col6: lca_rank (index 6)
            taxid = cols[2]
            lca_rank = cols[6]
            
            if taxid not in aggregates:
                aggregates[taxid] = {'total': 0, 'species': 0, 'genus': 0}
            
            aggregates[taxid]['total'] += 1
            if lca_rank == 'species':
                aggregates[taxid]['species'] += 1
            elif lca_rank == 'genus':
                aggregates[taxid]['genus'] += 1

    # 4 & 5. Read coverage.tsv and assign columns/flag
    with open(args.coverage, 'r') as fin, open(args.output, 'w') as fout:
        header_line = fin.readline().rstrip("\n")
        header = header_line.split("\t")
        # Add new headers
        new_header = header + ["Total", "Species", "Genus", "Flag"]
        fout.write("\t".join(new_header) + "\n")
        
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split('\t')
            # Taxon, Status, Taxid, Color
            # Taxid is cols[2]
            if len(cols) < 3:
                fout.write(line + "\n")
                continue
                
            status = cols[1]
            taxid = cols[2]
            
            if status == "OK":
                if taxid in aggregates:
                    a = aggregates[taxid]
                    total = str(a['total'])
                    species = str(a['species'])
                    genus = str(a['genus'])
                    flag = "-"
                else:
                    total = species = genus = "0"
                    flag = "member_only"
            else:
                total = species = genus = "0"
                flag = "-"
            
            out_cols = cols + [total, species, genus, flag]
            fout.write("\t".join(out_cols) + "\n")

if __name__ == "__main__":
    main()
