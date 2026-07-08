#!/usr/bin/env python
import sys

IUPAC = {
    frozenset('A'): 'A', frozenset('C'): 'C', frozenset('G'): 'G', frozenset('T'): 'T',
    frozenset('AG'): 'R', frozenset('CT'): 'Y', frozenset('CG'): 'S', frozenset('AT'): 'W', 
    frozenset('GT'): 'K', frozenset('AC'): 'M', frozenset('CGT'): 'B', frozenset('AGT'): 'D', 
    frozenset('ACT'): 'H', frozenset('ACG'): 'V', frozenset('ACGT'): 'N'
}

def collapse(seqs):
    # Combine the letters at each position into a single IUPAC degenerate code
    return "".join(IUPAC.get(frozenset(seq[i].upper() for seq in seqs), 'N') for i in range(len(seqs[0])))

def process_sample_sheet(input_file, output_file):
    # Read line by line
    with open(input_file, 'r') as f:
        lines = [line.strip().split('\t') for line in f if line.strip()]
        
    header, data = lines[0], lines[1:]
    
    # Group by primer name
    groups = {}
    for row in data:
        groups.setdefault(row[0], []).append(row)
        
    out = [header]
    for primer, rows in groups.items():
        # If duplicated, check if all forward lengths and all reverse lengths are the same
        if len(rows) > 1 and all(len(r[1]) == len(rows[0][1]) for r in rows) and all(len(r[2]) == len(rows[0][2]) for r in rows):
            # Lengths match: collapse and replace
            fwd = collapse([r[1] for r in rows])
            rev = collapse([r[2] for r in rows])
            out.append([primer, fwd, rev, rows[0][3], rows[0][4]])
        else:
            # Lengths don't match or not a duplicate: do nothing, just keep original
            out.extend(rows)
            
    # Write everything back
    with open(output_file, 'w') as f:
        f.write('\n'.join('\t'.join(row) for row in out) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: process_sample_sheet.py <input.tsv> <output.tsv>")
        sys.exit(1)
    process_sample_sheet(sys.argv[1], sys.argv[2])
