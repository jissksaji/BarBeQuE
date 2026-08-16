#!/usr/bin/env python3
"""Collapse a FASTA of forward/reverse primer variants into one consensus
forward primer and one consensus reverse primer.

Direction is inferred from an `fwd`/`forward`/`rev`/`reverse` token in each
header (case-insensitive), matching the naming convention already used for
primer FASTAs in this ecosystem (e.g. MA_FWD, POL_REV, Fwd 1.1). All
fwd-labelled sequences are collapsed together into one IUPAC-degenerate
consensus, and all rev-labelled sequences are collapsed together into
another - regardless of any name-prefix grouping (e.g. MA vs POL) - since
obipcr/cutadapt each only accept one forward and one reverse primer string.

Collapsing (building an IUPAC consensus) only makes sense when 2+ variants
are given for the same direction - many real primer FASTAs are just a single,
already-specific fwd/rev pair with no fwd/rev token at all (e.g. rbcL1/rbcLB,
CLP1a/CLP2). For that unambiguous case - exactly 2 records, no direction
tokens anywhere - the 1st record is treated as fwd and the 2nd as rev, with
nothing to collapse. Anything else untagged (1 record, or 3+ with no tokens,
or a mix of tagged and untagged) is a hard error rather than a guess.
"""

import argparse
import re
import sys

# Same consensus table as bin/process_sample_sheet.py
IUPAC = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M", frozenset("CGT"): "B", frozenset("AGT"): "D",
    frozenset("ACT"): "H", frozenset("ACG"): "V", frozenset("ACGT"): "N",
}

DIRECTION_RE = re.compile(r"(?i)(fwd|forward|rev|reverse)")


def read_fasta(path):
    """Yield (header, sequence) pairs, joining multi-line sequences."""
    header = None
    seq_parts = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def classify(header):
    match = DIRECTION_RE.search(header)
    if not match:
        return None
    return "fwd" if match.group(1).lower() in ("fwd", "forward") else "rev"


EXPAND = {v: k for k, v in IUPAC.items()}


def collapse(seqs):
    """Combine the letters at each position into a single IUPAC degenerate code."""
    return "".join(
        IUPAC.get(frozenset().union(*(EXPAND.get(seq[i].upper(), frozenset("ACGT")) for seq in seqs)), "N")
        for i in range(len(seqs[0]))
    )


def resolve(seqs, direction, fasta_path):
    if len({len(s) for s in seqs}) != 1:
        sys.exit(
            f"{direction} primer variants in {fasta_path} have mismatched lengths "
            f"({[len(s) for s in seqs]}) - cannot collapse into one consensus sequence."
        )
    return seqs[0] if len(seqs) == 1 else collapse(seqs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="FASTA with fwd/rev-labelled headers")
    parser.add_argument("--prefix", required=True, help="Prefix for the two output record IDs")
    parser.add_argument("--out", required=True, help="Output FASTA path")
    args = parser.parse_args()

    records = list(read_fasta(args.fasta))
    classified = [(header, seq, classify(header)) for header, seq in records]
    untagged = [header for header, _, direction in classified if direction is None]

    groups = {"fwd": [], "rev": []}
    if not untagged:
        for _header, seq, direction in classified:
            groups[direction].append(seq)
    elif len(untagged) == len(records) and len(records) == 2:
        # No fwd/rev tokens anywhere, and exactly one record per direction - unambiguous,
        # nothing to collapse.
        groups["fwd"].append(records[0][1])
        groups["rev"].append(records[1][1])
    else:
        sys.exit(
            f"{args.fasta} has header(s) with no fwd/forward/rev/reverse token "
            f"({untagged}) and isn't a plain 2-record fwd/rev pair - cannot tell "
            "which direction each sequence belongs to."
        )

    if not groups["fwd"] or not groups["rev"]:
        sys.exit(f"{args.fasta} is missing a fwd or rev sequence - found {len(groups['fwd'])} fwd, {len(groups['rev'])} rev.")

    fwd = resolve(groups["fwd"], "fwd", args.fasta)
    rev = resolve(groups["rev"], "rev", args.fasta)

    with open(args.out, "w") as handle:
        handle.write(f">{args.prefix}_fwd\n{fwd}\n")
        handle.write(f">{args.prefix}_rev\n{rev}\n")


if __name__ == "__main__":
    main()
