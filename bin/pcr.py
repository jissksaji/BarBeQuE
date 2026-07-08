#!/usr/bin/env python3

import argparse
import sys
from bisect import bisect_left, bisect_right
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from numba import njit


MASK_TABLE = np.zeros(256, dtype=np.uint8)

for c in b"Aa":
    MASK_TABLE[c] = 1
for c in b"Cc":
    MASK_TABLE[c] = 2
for c in b"Gg":
    MASK_TABLE[c] = 4
for c in b"TtUu":
    MASK_TABLE[c] = 8
for c in b"Rr":
    MASK_TABLE[c] = 1 | 4
for c in b"Yy":
    MASK_TABLE[c] = 2 | 8
for c in b"Ss":
    MASK_TABLE[c] = 2 | 4
for c in b"Ww":
    MASK_TABLE[c] = 1 | 8
for c in b"Kk":
    MASK_TABLE[c] = 4 | 8
for c in b"Mm":
    MASK_TABLE[c] = 1 | 2
for c in b"Bb":
    MASK_TABLE[c] = 2 | 4 | 8
for c in b"Dd":
    MASK_TABLE[c] = 1 | 4 | 8
for c in b"Hh":
    MASK_TABLE[c] = 1 | 2 | 8
for c in b"Vv":
    MASK_TABLE[c] = 1 | 2 | 4
for c in b"Nn":
    MASK_TABLE[c] = 0


COMP = str.maketrans(
    "ACGTRYSWKMBDHVNUacgtryswkmbdhvnu",
    "TGCAYRSWMKVHDBNAtgcayrswmkvhdbna",
)


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def read_fasta(path):
    name = None
    parts = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()

                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

        if name is not None:
            yield name, "".join(parts).upper()


def to_mask_array(seq: str) -> np.ndarray:
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    return MASK_TABLE[arr]


@njit(cache=True, nogil=True)
def find_hits_numba(seq_masks, primer_masks, max_mismatches):
    n = seq_masks.shape[0]
    k = primer_masks.shape[0]

    if n < k:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int16)

    hit_count = 0

    for i in range(n - k + 1):
        mismatches = 0

        for j in range(k):
            if seq_masks[i + j] & primer_masks[j] == 0:
                mismatches += 1

                if mismatches > max_mismatches:
                    break

        if mismatches <= max_mismatches:
            hit_count += 1

    starts = np.empty(hit_count, dtype=np.int64)
    errors = np.empty(hit_count, dtype=np.int16)

    idx = 0

    for i in range(n - k + 1):
        mismatches = 0

        for j in range(k):
            if seq_masks[i + j] & primer_masks[j] == 0:
                mismatches += 1

                if mismatches > max_mismatches:
                    break

        if mismatches <= max_mismatches:
            starts[idx] = i
            errors[idx] = mismatches
            idx += 1

    return starts, errors


def write_fasta(header, seq, width=0, stream=None):
    if stream is None:
        stream = sys.stdout
    print(f">{header}", file=stream)

    if width == 0:
        print(seq, file=stream)
    else:
        for i in range(0, len(seq), width):
            print(seq[i:i + width], file=stream)


def extract_from_sequence(
    seq_id,
    seq,
    fwd,
    rev,
    mismatches,
    min_len,
    max_len,
    include_primers=False,
):
    """Find all amplicons in one sequence. Returns a list of (header, amplicon)
    tuples instead of printing directly, so this can safely be called from
    multiple threads at once -- the caller writes the results out."""
    rev_rc = reverse_complement(rev)

    fwd_masks = to_mask_array(fwd)
    rev_rc_masks = to_mask_array(rev_rc)

    fwd_len = len(fwd)
    rev_len = len(rev_rc)

    results = []

    for direction, strand_seq in [
        ("forward", seq),
        ("reverse", reverse_complement(seq)),
    ]:
        strand_masks = to_mask_array(strand_seq)

        f_starts, f_errors = find_hits_numba(
            strand_masks,
            fwd_masks,
            mismatches,
        )

        r_starts, r_errors = find_hits_numba(
            strand_masks,
            rev_rc_masks,
            mismatches,
        )

        if len(f_starts) == 0 or len(r_starts) == 0:
            continue

        r_starts_list = r_starts.tolist()

        for f_idx, f_start in enumerate(f_starts):
            f_start = int(f_start)
            f_end = f_start + fwd_len

            min_r_start = f_end + min_len
            max_r_start = f_end + max_len

            lo = bisect_left(r_starts_list, min_r_start)
            hi = bisect_right(r_starts_list, max_r_start)

            # take only the NEAREST valid end -- avoids reporting a spurious
            # "merged" amplicon that spans across two separate real occurrences
            if lo < hi:
                r_idx = lo
                r_start = int(r_starts[r_idx])
                r_end = r_start + rev_len

                if include_primers:
                    amp_start = f_start
                    amp_end = r_end
                else:
                    amp_start = f_end
                    amp_end = r_start

                amplicon = strand_seq[amp_start:amp_end].upper()

                header = (
                    f"{seq_id}_sub[{amp_start + 1}..{amp_end}]"
                    f" direction={direction}"
                    f" forward_error={int(f_errors[f_idx])}"
                    f" reverse_error={int(r_errors[r_idx])}"
                    f" forward_match={strand_seq[f_start:f_end]}"
                    f" reverse_match={strand_seq[r_start:r_end]}"
                )

                results.append((header, amplicon))

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Fast sliding-window in-silico PCR using Numba. No indels. Integer mismatches only."
    )

    parser.add_argument("-i", "--input", required=True, help="Input FASTA")
    parser.add_argument("-f", "--forward", required=True, help="Forward primer 5'->3'")
    parser.add_argument("-r", "--reverse", required=True, help="Reverse primer 5'->3'")
    parser.add_argument("-m", "--mismatches", type=int, default=0,
                        help="Allowed mismatches per primer")
    parser.add_argument("--min-length", type=int, default=0,
                        help="Minimum internal amplicon length, primers excluded")
    parser.add_argument("--max-length", type=int, required=True,
                        help="Maximum internal amplicon length, primers excluded")
    parser.add_argument("--include-primers", action="store_true",
                        help="Output amplicons including primer-binding sites")
    parser.add_argument("-w", "--width", type=int, default=0,
                        help="FASTA line width. 0 = one-line sequence")
    parser.add_argument("-o", "--output", default=None,
                        help="Output FASTA file (default: stdout)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of worker threads (default: 1)")

    args = parser.parse_args()

    fwd = args.forward.upper()
    rev = args.reverse.upper()

    out_stream = open(args.output, 'w') if args.output else sys.stdout

    def process_one(record):
        seq_id, seq = record
        return extract_from_sequence(
            seq_id=seq_id,
            seq=seq,
            fwd=fwd,
            rev=rev,
            mismatches=args.mismatches,
            min_len=args.min_length,
            max_len=args.max_length,
            include_primers=args.include_primers,
        )

    try:
        if args.threads <= 1:
            # single-threaded path -- no pool overhead
            for record in read_fasta(args.input):
                for header, amplicon in process_one(record):
                    write_fasta(header, amplicon, width=args.width, stream=out_stream)
        else:
            # Bounded in-flight queue: only read ahead a few records at a time,
            # rather than loading the entire (potentially huge) FASTA into
            # memory before any work starts.
            max_in_flight = args.threads * 4
            with ThreadPoolExecutor(max_workers=args.threads) as pool:
                record_iter = read_fasta(args.input)
                in_flight = {}  # future -> True, just used as an ordered set

                def fill_queue():
                    for record in record_iter:
                        in_flight[pool.submit(process_one, record)] = True
                        if len(in_flight) >= max_in_flight:
                            return True
                    return False

                more_to_read = fill_queue()
                while in_flight:
                    finished = next(as_completed(in_flight))
                    del in_flight[finished]
                    for header, amplicon in finished.result():
                        write_fasta(header, amplicon, width=args.width, stream=out_stream)
                    if more_to_read:
                        more_to_read = fill_queue()
    finally:
        if args.output:
            out_stream.close()


if __name__ == "__main__":
    main()