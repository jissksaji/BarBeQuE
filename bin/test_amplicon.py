#!/usr/bin/env python3
"""
find_all_amplicons.py

Recovers ALL primer-pair occurrences in EVERY orientation from a single
input sequence, in ONE pass per record -- no iteration, no masking needed.

Adapted from AmpliconHunter's core algorithm (Hyperscan-based all-matches
scan + exhaustive F/R pairing), reimplemented here with Python's `regex`
module for fuzzy/degenerate matching, so it runs against a normal
multi-record reference FASTA (your actual database format) rather than
AmpliconHunter's one-genome-per-file convention.

ALGORITHM
---------
1. For each sequence, find every occurrence of all four patterns
   (forward primer, its reverse complement, reverse primer, its reverse
   complement) using substitution-only fuzzy matching (mirrors --no-indels).
2. Pair every F/R "start" match with every compatible "end" match
   downstream within [min_length, max_length] -- exactly AmpliconHunter's
   nested-loop pairing logic. This naturally finds:
     - multiple amplicons in the same record (no masking required)
     - amplicons in both the given orientation AND on the opposite strand
       (an RF/RR-orientation hit is reverse-complemented back to sense
       before output, same normalization AmpliconHunter does)
3. Optional 3'-clamp: reject a match if its 3'-most N bases (the end most
   sensitive to real PCR extension) contain any substitution, even if the
   overall mismatch budget would allow it.

USAGE
-----
  python3 find_all_amplicons.py \\
      --input refs.fasta --output amplicons.fasta \\
      --fwd ATGCGATACTTGGTGTGAAT --rev TCCTCCGCTTATTGATATGC \\
      --mismatches 2 --clamp 0 --min-length 1 --max-length 1000

Output records are named '<original_id>_amp1', '_amp2', ... when a record
yields more than one amplicon; a single-hit record keeps its original name.
"""

import argparse
import sys
import regex

IUPAC_COMP = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D', 'N': 'N',
}
IUPAC_CLASS = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[CG]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}


def revcomp(seq):
    return ''.join(IUPAC_COMP[b] for b in reversed(seq.upper()))


def to_fuzzy_pattern(primer, mismatches):
    """Build a substitution-only fuzzy regex for a (possibly degenerate) primer."""
    body = ''.join(IUPAC_CLASS[b] for b in primer.upper())
    if mismatches == 0:
        return regex.compile(body)
    return regex.compile(f'(?:{body}){{s<={mismatches}}}')


def find_all_matches(seq, pattern, label):
    """Return list of (label, start, end, n_subs) for every non-overlapping match."""
    results = []
    for m in regex.finditer(pattern, seq, overlapped=False):
        n_subs = m.fuzzy_counts[0] if m.fuzzy_counts else 0
        results.append((label, m.start(), m.end(), n_subs))
    return results


def passes_clamp(seq, start, end, primer_literal, label, clamp):
    """Check the 3'-most `clamp` bases of the matched primer for exact agreement.
    For F/R (primer read left-to-right as given), the 3' end is the right edge.
    For F_rc/R_rc (primer matched on its reverse complement), the 3' end of the
    ORIGINAL primer corresponds to the LEFT edge of this match."""
    if clamp <= 0:
        return True
    matched = seq[start:end]
    if label in ('F', 'R'):
        primer_tail = primer_literal[-clamp:]
        seq_tail = matched[-clamp:]
    else:  # F_rc, R_rc -- the primer's 3' end is at the start of the rc match
        primer_tail = revcomp(primer_literal)[:clamp]
        seq_tail = matched[:clamp]
    for p_base, s_base in zip(primer_tail.upper(), seq_tail.upper()):
        allowed = IUPAC_CLASS[p_base].strip('[]')
        if s_base not in allowed:
            return False
    return True


def read_fasta(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, ''.join(chunks)


def find_amplicons_in_sequence(seq, fwd, rev, mismatches, clamp, min_len, max_len, include_offtarget):
    fwd_rc = revcomp(fwd)
    rev_rc = revcomp(rev)

    patterns = {
        'F': to_fuzzy_pattern(fwd, mismatches),
        'F_rc': to_fuzzy_pattern(fwd_rc, mismatches),
        'R': to_fuzzy_pattern(rev, mismatches),
        'R_rc': to_fuzzy_pattern(rev_rc, mismatches),
    }
    primer_for_label = {'F': fwd, 'F_rc': fwd, 'R': rev, 'R_rc': rev}

    all_matches = []
    for label, pattern in patterns.items():
        all_matches.extend(find_all_matches(seq, pattern, label))
    all_matches.sort(key=lambda m: m[1])  # sort by start position

    amplicons = []  # list of (interior_seq, orientation)
    used_end_indices = set()

    for i, (label_i, start_i, end_i, _) in enumerate(all_matches):
        if label_i not in ('F', 'R'):
            continue
        for j in range(i + 1, len(all_matches)):
            label_j, start_j, end_j, _ = all_matches[j]
            length = end_j - start_i
            if length > max_len:
                break
            if label_j in ('F', 'R'):
                break  # hit another start before finding a valid end -- stop this branch
            if j in used_end_indices:
                continue
            if length < min_len:
                continue
            orientation = label_i[0] + label_j[0]  # e.g. 'F'+'R' -> 'FR', 'R'+'F' -> 'RF'
            if not include_offtarget and orientation not in ('FR', 'RF'):
                continue

            if not passes_clamp(seq, start_i, end_i, primer_for_label[label_i], label_i, clamp):
                continue
            if not passes_clamp(seq, start_j, end_j, primer_for_label[label_j], label_j, clamp):
                continue

            interior = seq[end_i:start_j]
            if orientation[0] == 'R':
                interior = revcomp(interior)

            amplicons.append((interior, orientation))
            used_end_indices.add(j)
            break  # this start is consumed; move to the next start match

    return amplicons


def main():
    ap = argparse.ArgumentParser(description="Find all primer-pair amplicons in every orientation, single pass per record.")
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--fwd', required=True)
    ap.add_argument('--rev', required=True)
    ap.add_argument('--mismatches', type=int, required=True)
    ap.add_argument('--clamp', type=int, default=0, help="3'-most N bases of each primer must match exactly (default 0 = off)")
    ap.add_argument('--min-length', type=int, default=1)
    ap.add_argument('--max-length', type=int, default=10000)
    ap.add_argument('--include-offtarget', action='store_true', help="also keep FF/RR orientation hits")
    args = ap.parse_args()

    fwd = args.fwd.upper()
    rev = args.rev.upper()

    n_in, n_out, n_multi = 0, 0, 0
    with open(args.output, 'w') as out_f:
        for name, seq in read_fasta(args.input):
            n_in += 1
            amplicons = find_amplicons_in_sequence(
                seq.upper(), fwd, rev, args.mismatches, args.clamp,
                args.min_length, args.max_length, args.include_offtarget,
            )
            if len(amplicons) > 1:
                n_multi += 1
            for idx, (interior, orientation) in enumerate(amplicons, 1):
                out_name = name if len(amplicons) == 1 else f'{name}_amp{idx}'
                out_f.write(f'>{out_name} orientation={orientation}\n{interior}\n')
                n_out += 1

    print(f"Input records:          {n_in}", file=sys.stderr)
    print(f"Output amplicons:       {n_out}", file=sys.stderr)
    print(f"Records with >1 hit:    {n_multi}", file=sys.stderr)


if __name__ == '__main__':
    main()