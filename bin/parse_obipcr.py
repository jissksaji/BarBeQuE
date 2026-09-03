#!/usr/bin/env python3
import sys
import json
import re
import argparse

# Expected GC contribution of each IUPAC nucleotide.  Fractional values make
# GC% meaningful for collapsed/degenerate primers while remaining identical
# to the usual calculation for unambiguous A/C/G/T sequences.
IUPAC_GC_FRACTION = {
    'A': 0.0, 'C': 1.0, 'G': 1.0, 'T': 0.0, 'U': 0.0,
    'R': 0.5, 'Y': 0.5, 'S': 1.0, 'W': 0.0, 'K': 0.5, 'M': 0.5,
    'B': 2.0 / 3.0, 'D': 1.0 / 3.0, 'H': 1.0 / 3.0,
    'V': 2.0 / 3.0, 'N': 0.5,
}

def parse_fasta_header(header_line):
    """Parses a single OBI ecoPCR fasta header."""
    header_line = header_line.strip()
    if not header_line.startswith(">"):
        return None
    
    parts = header_line.split(" ", 1)
    if len(parts) != 2:
        return None
        
    seq_info = parts[0][1:]
    json_str = parts[1]
    
    match = re.match(r"(.+)_sub\[(\d+)\.\.(\d+)\]", seq_info)
    if not match:
        return None
        
    seq_id = match.group(1)
    start_pos = int(match.group(2))
    end_pos = int(match.group(3))
    
    try:
        metadata = json.loads(json_str)
    except json.JSONDecodeError:
        return None
        
    return {
        "seq_id": seq_id,
        "start_pos": start_pos,
        "end_pos": end_pos,
        "metadata": metadata
    }

def calculate_gc(sequence):
    """Calculate expected GC%, including fractional IUPAC ambiguity."""
    if not sequence: return "0.0"
    seq_upper = sequence.upper()
    gc_count = sum(IUPAC_GC_FRACTION.get(base, 0.0) for base in seq_upper)
    return str(round((gc_count / len(sequence)) * 100, 2))

IUPAC_CODES = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'U', 'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}
}

def is_iupac_match(primer_base, match_base):
    """Return whether a primer pattern base matches a reference base.

    IUPAC ambiguity is directional here: ambiguity codes in the primer expand
    the bases it can bind, while an ambiguous base in the reference sequence
    is counted as an OBI-PCR error because its actual nucleotide is unknown.
    """
    primer_bases = IUPAC_CODES.get(
        primer_base.upper(),
        {primer_base.upper()},
    )
    match_base = match_base.upper()
    return match_base in {"A", "C", "G", "T"} and match_base in primer_bases


def mismatch_positions(primer, match_seq):
    """Return 1-based mismatch positions from the primer's 5' and 3' ends."""
    primer = primer.upper()
    match_seq = match_seq.upper()
    length = min(len(primer), len(match_seq))
    indices = [
        i
        for i in range(length)
        if not is_iupac_match(primer[i], match_seq[i])
    ]
    return (
        [i + 1 for i in indices],
        [length - i for i in indices],
    )

def process_obipcr(input_file, output_file):
    records = []

    with open(input_file, 'r') as infile:
        current_header = None
        current_seq = []

        def flush_record(header_line, seq):
            parsed = parse_fasta_header(header_line)
            if not parsed:
                return

            seq_id = parsed["seq_id"]
            start_pos = parsed["start_pos"]
            end_pos = parsed["end_pos"]
            amplicon_length = end_pos - start_pos + 1
            meta = parsed["metadata"]

            fw_errors = meta.get("forward_error", 0)
            rv_errors = meta.get("reverse_error", 0)
            # obipcr echoes the primer exactly as it was given, so it still carries the
            # '#' 3'-clamp markers added by --obipcr_fixed_3prime. Strip them, otherwise
            # they shift every position against the match and corrupt the primer metrics.
            fw_primer = meta.get("forward_primer", "").replace("#", "")
            fw_match = meta.get("forward_match", "")
            rv_primer = meta.get("reverse_primer", "").replace("#", "")
            rv_match = meta.get("reverse_match", "")

            fw_mm_positions, fw_mm_3prime_positions = mismatch_positions(fw_primer, fw_match)
            rv_mm_positions, rv_mm_3prime_positions = mismatch_positions(rv_primer, rv_match)
            fw_mm_primer = ",".join(map(str, fw_mm_positions))
            rv_mm_primer = ",".join(map(str, rv_mm_positions))
            fw_mm_3prime = ",".join(map(str, fw_mm_3prime_positions))
            rv_mm_3prime = ",".join(map(str, rv_mm_3prime_positions))

            row = [
                seq_id,
                str(amplicon_length),
                fw_primer,
                fw_match,
                str(fw_errors),
                fw_mm_primer if fw_mm_primer else "None",
                fw_mm_3prime if fw_mm_3prime else "None",
                calculate_gc(fw_primer),
                calculate_gc(fw_match),

                rv_primer,
                rv_match,
                str(rv_errors),
                rv_mm_primer if rv_mm_primer else "None",
                rv_mm_3prime if rv_mm_3prime else "None",
                calculate_gc(rv_primer),
                calculate_gc(rv_match),
                calculate_gc(seq),
                seq
            ]
            records.append({"seq_id": seq_id, "amplicon_length": amplicon_length, "row": row})

        for line in infile:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if current_header:
                    flush_record(current_header, "".join(current_seq))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        if current_header:
            flush_record(current_header, "".join(current_seq))

    with open(output_file, 'w') as outfile:
        headers = [
            "Sequence_ID",
            "Amplicon_Length",
            "Forward_Primer",
            "Forward_Match",
            "Forward_Errors",
            "Forward_Mismatch_Positions_Primer",
            "Forward_Mismatch_From_3Prime",
            "Forward_Primer_GC",
            "Forward_Match_GC",
            "Reverse_Primer",
            "Reverse_Match",
            "Reverse_Errors",
            "Reverse_Mismatch_Positions_Primer",
            "Reverse_Mismatch_From_3Prime",
            "Reverse_Primer_GC",
            "Reverse_Match_GC",
            "Amplicon_GC",
            "Amplicon_Sequence"
        ]
        outfile.write("\t".join(headers) + "\n")

        for r in records:
            outfile.write("\t".join(r["row"]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse OBI ecoPCR output to TSV.")
    parser.add_argument("input", help="Input fasta-like file from ecoPCR")
    parser.add_argument("output", help="Output TSV file")
    args = parser.parse_args()
    
    process_obipcr(args.input, args.output)
    print(f"Extraction complete. Results saved to {args.output}")
