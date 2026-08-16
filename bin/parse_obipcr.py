#!/usr/bin/env python3
import sys
import json
import re
import argparse
import itertools
from collections import defaultdict

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

def calculate_tm(sequence):
    """
    Estimates primer:template melting temperature (Tm) in degrees C.
    Uses the Wallace rule for short oligos (<14 nt) and the salt-adjusted
    GC formula otherwise. Computed on the bound (match) sequence, so it
    reflects the duplex that actually forms rather than the degenerate primer.
    """
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return "0.0"
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')
    if n < 14:
        tm = 2.0 * at + 4.0 * gc
    else:
        tm = 64.9 + 41.0 * (gc - 16.4) / n
    return str(round(tm, 1))

def gc_clamp_count(sequence):
    """Returns number of G/C in the last 5 bases of 3' end."""
    if len(sequence) < 5: return "0"
    end = sequence[-5:].upper()
    return str(end.count('G') + end.count('C'))

def max_mononucleotide_run(sequence):
    """Finds the longest run of a single nucleotide."""
    if not sequence: return "0"
    return str(max((sum(1 for _ in group) for key, group in itertools.groupby(sequence.upper())), default=0))

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, 'N') for base in reversed(seq.upper()))

def three_prime_dimer_length(seq1, seq2):
    """
    Longest complementary run anchored at the 3' terminal base of seq1 against seq2,
    modeling a primer-dimer that polymerase could extend from (the biologically
    relevant dimer risk, unlike a plain longest-common-substring-anywhere check).
    """
    if not seq1 or not seq2: return "0"
    seq1 = seq1.upper()
    rc2 = reverse_complement(seq2)
    best = 0
    for offset in range(-(len(rc2) - 1), len(seq1)):
        i = len(seq1) - 1
        j = i - offset
        run = 0
        while i >= 0 and 0 <= j < len(rc2) and seq1[i] == rc2[j]:
            run += 1
            i -= 1
            j -= 1
        best = max(best, run)
    return str(best)

def max_hairpin_stem(sequence, min_loop=3):
    """
    Longest self-complementary stem forming a fold-back hairpin, requiring an
    unpaired loop of at least min_loop bases so the fold is physically possible
    (unlike testing full reverse-complement self-similarity with no loop constraint).
    """
    if not sequence: return "0"
    seq = sequence.upper()
    n = len(seq)
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    best = 0
    for i in range(n):
        for j in range(i + min_loop + 1, n):
            if complement.get(seq[i]) != seq[j]:
                continue
            stem = 0
            li, rj = i, j
            while li >= 0 and rj < n and complement.get(seq[li]) == seq[rj]:
                stem += 1
                li -= 1
                rj += 1
            best = max(best, stem)
    return str(best)

MISMATCH_SEVERITY = {
    frozenset({'G', 'T'}): 1.0,
    frozenset({'A', 'C'}): 1.5,
    frozenset({'C', 'T'}): 1.5,
    frozenset({'A', 'G'}): 2.0,
    frozenset({'A', 'T'}): 2.0,
    frozenset({'C', 'G'}): 2.0,
}

IUPAC_CODES = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'U', 'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}
}

def is_iupac_match(b1, b2):
    b1_bases = IUPAC_CODES.get(b1.upper(), {b1.upper()})
    b2_bases = IUPAC_CODES.get(b2.upper(), {b2.upper()})
    return bool(b1_bases & b2_bases)

def find_mismatches(primer, match_seq, start_pos, is_reverse=False):
    """
    Finds mismatches and returns coordinates plus a severity-weighted score:
    - 1-based index from 5' end
    - absolute genomic position
    - 1-based index from 3' end
    - severity score: base-pair destabilization weight, amplified for
      mismatches near the 3' end since those are far more likely to block
      polymerase extension than 5' mismatches.
    """
    primer = primer.upper()
    match_seq = match_seq.upper()

    primer_mismatches = []
    genome_mismatches = []
    three_prime_mismatches = []
    severity_total = 0.0

    length = min(len(primer), len(match_seq))
    for i in range(length):
        if not is_iupac_match(primer[i], match_seq[i]):
            dist_from_3prime = length - i
            primer_mismatches.append(str(i + 1))
            three_prime_mismatches.append(str(dist_from_3prime))

            if not is_reverse:
                genome_pos = start_pos + i
            else:
                genome_pos = start_pos - i

            genome_mismatches.append(str(genome_pos))

            base_severity = MISMATCH_SEVERITY.get(frozenset({primer[i], match_seq[i]}), 2.0)
            if dist_from_3prime <= 3:
                position_weight = 3.0
            elif dist_from_3prime <= 6:
                position_weight = 1.5
            else:
                position_weight = 1.0
            severity_total += base_severity * position_weight

    return (",".join(primer_mismatches), ",".join(genome_mismatches),
            ",".join(three_prime_mismatches), str(round(severity_total, 2)))

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

            direction = meta.get("direction", "")
            fw_errors = meta.get("forward_error", 0)
            rv_errors = meta.get("reverse_error", 0)
            # obipcr echoes the primer exactly as it was given, so it still carries the
            # '#' 3'-clamp markers added by --obipcr_fixed_3prime. Strip them, otherwise
            # they shift every position against the match and corrupt the primer metrics.
            fw_primer = meta.get("forward_primer", "").replace("#", "")
            fw_match = meta.get("forward_match", "")
            rv_primer = meta.get("reverse_primer", "").replace("#", "")
            rv_match = meta.get("reverse_match", "")

            fw_mm_primer, fw_mm_genome, fw_mm_3prime, fw_mm_severity = find_mismatches(fw_primer, fw_match, start_pos, is_reverse=False)
            rv_mm_primer, rv_mm_genome, rv_mm_3prime, rv_mm_severity = find_mismatches(rv_primer, rv_match, end_pos, is_reverse=True)

            amplicon_gc = calculate_gc(seq)
            hetero_dimer = max(
                int(three_prime_dimer_length(fw_primer, rv_primer)),
                int(three_prime_dimer_length(rv_primer, fw_primer))
            )

            fw_tm = calculate_tm(fw_match)
            rv_tm = calculate_tm(rv_match)
            tm_diff = str(round(abs(float(fw_tm) - float(rv_tm)), 1))

            row = [
                seq_id,
                str(start_pos),
                str(end_pos),
                str(amplicon_length),
                direction,

                fw_primer,
                fw_match,
                str(fw_errors),
                fw_mm_primer if fw_mm_primer else "None",
                fw_mm_3prime if fw_mm_3prime else "None",
                fw_mm_genome if fw_mm_genome else "None",
                fw_mm_severity,
                calculate_gc(fw_primer),
                calculate_gc(fw_match),
                fw_tm,
                gc_clamp_count(fw_primer),
                max_mononucleotide_run(fw_primer),
                three_prime_dimer_length(fw_primer, fw_primer),
                max_hairpin_stem(fw_primer),

                rv_primer,
                rv_match,
                str(rv_errors),
                rv_mm_primer if rv_mm_primer else "None",
                rv_mm_3prime if rv_mm_3prime else "None",
                rv_mm_genome if rv_mm_genome else "None",
                rv_mm_severity,
                calculate_gc(rv_primer),
                calculate_gc(rv_match),
                rv_tm,
                gc_clamp_count(rv_primer),
                max_mononucleotide_run(rv_primer),
                three_prime_dimer_length(rv_primer, rv_primer),
                max_hairpin_stem(rv_primer),

                str(hetero_dimer),
                tm_diff,
                amplicon_gc,
                max_mononucleotide_run(seq),
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

    hits_by_seq_id = defaultdict(list)
    for r in records:
        hits_by_seq_id[r["seq_id"]].append(r["amplicon_length"])

    with open(output_file, 'w') as outfile:
        headers = [
            "Sequence_ID",
            "Forward_Binding_Start",
            "Reverse_Binding_End",
            "Amplicon_Length",
            "Direction",
            "Forward_Primer",
            "Forward_Match",
            "Forward_Errors",
            "Forward_Mismatch_Positions_Primer",
            "Forward_Mismatch_From_3Prime",
            "Forward_Mismatch_Positions_Genome",
            "Forward_Mismatch_Severity",
            "Forward_Primer_GC",
            "Forward_Match_GC",
            "Forward_Tm",
            "Forward_GC_Clamp",
            "Forward_Max_Run",
            "Forward_Self_Dimer_3prime",
            "Forward_Hairpin_Stem",
            "Reverse_Primer",
            "Reverse_Match",
            "Reverse_Errors",
            "Reverse_Mismatch_Positions_Primer",
            "Reverse_Mismatch_From_3Prime",
            "Reverse_Mismatch_Positions_Genome",
            "Reverse_Mismatch_Severity",
            "Reverse_Primer_GC",
            "Reverse_Match_GC",
            "Reverse_Tm",
            "Reverse_GC_Clamp",
            "Reverse_Max_Run",
            "Reverse_Self_Dimer_3prime",
            "Reverse_Hairpin_Stem",
            "Hetero_Dimer_3prime_Len",
            "Tm_Difference",
            "Amplicon_GC",
            "Amplicon_Max_Run",
            "Amplicon_Sequence",
            "Hits_On_Sequence",
            "Amplicon_Length_Spread"
        ]
        outfile.write("\t".join(headers) + "\n")

        for r in records:
            lengths = hits_by_seq_id[r["seq_id"]]
            hit_count = len(lengths)
            spread = max(lengths) - min(lengths)
            row = r["row"] + [str(hit_count), str(spread)]
            outfile.write("\t".join(row) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse OBI ecoPCR output to TSV.")
    parser.add_argument("input", help="Input fasta-like file from ecoPCR")
    parser.add_argument("output", help="Output TSV file")
    args = parser.parse_args()
    
    process_obipcr(args.input, args.output)
    print(f"Extraction complete. Results saved to {args.output}")
