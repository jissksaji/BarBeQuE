#!/usr/bin/env python3
"""Turn primer FASTA input into the TSV samplesheet the pipeline already reads.

`--input` is either a directory of primer FASTAs or a single primer FASTA.
Every file is parsed into one or more `primer / fwd / rev / min / max` rows, so
all three input forms (samplesheet, single FASTA, folder of FASTAs) converge on
the same shape before anything else runs.

Records are grouped by the *prefix* of their header - the text in front of an
`fwd`/`forward`/`rev`/`reverse` token. Anything after the token is a variant
label and is ignored, so `MA_FWD_1` and `MA_FWD_2` are both prefix `MA`. Each
prefix becomes its own primer set, which is what keeps two different markers in
one file (MA and POL) from being merged just because their primers happen to be
the same length.

Within a prefix, same-length variants collapse into one IUPAC-degenerate
primer, since obipcr takes a single forward and a single reverse string. When
the variants have different lengths they cannot be collapsed, so every forward
is instead paired with every reverse and each combination is benchmarked
separately - that is a warning, not an error.

A file is rejected only when it cannot yield a runnable pair at all: it holds no
FASTA records, a record has an empty sequence, a sequence contains a
non-nucleotide character, a prefix is missing one of its two directions, or a
record's direction cannot be determined. Every file is checked before anything
is written, so all problems are reported in one go.
"""

import argparse
import re
import sys
from pathlib import Path

# Same consensus table as bin/process_sample_sheet.py
IUPAC = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M", frozenset("CGT"): "B", frozenset("AGT"): "D",
    frozenset("ACT"): "H", frozenset("ACG"): "V", frozenset("ACGT"): "N",
}
EXPAND = {code: bases for bases, code in IUPAC.items()}
VALID_BASES = set(IUPAC.values()) | {"U"}

# The direction token must be delimited, so a prefix that merely contains
# "rev" or "fwd" (e.g. TREV_fwd) is not mistaken for the token itself.
DIRECTION_RE = re.compile(r"(?i)(?:^|[\s_.-])(fwd|forward|rev|reverse)(?:[\s_.-]|$)")
UNSAFE_NAME_RE = re.compile(r'[\s/\\:*?"<>|]')
FASTA_SUFFIXES = (".fa", ".fasta", ".fna")

TSV_HEADER = ["primer", "fwd", "rev", "min", "max"]


def read_fasta(path):
    """Yield (header, sequence) pairs, joining multi-line sequences."""
    header = None
    parts = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip()
                parts = []
            elif header is not None:
                parts.append(line)
        if header is not None:
            yield header, "".join(parts)


def split_header(header):
    """Return (prefix, direction), where direction is 'fwd', 'rev' or None."""
    match = DIRECTION_RE.search(header)
    if not match:
        return header, None
    direction = "fwd" if match.group(1).lower() in ("fwd", "forward") else "rev"
    return header[: match.start()].strip(" _.-"), direction


def collapse(seqs):
    """Combine the letters at each position into a single IUPAC degenerate code."""
    return "".join(
        IUPAC[frozenset().union(*(EXPAND[seq[i].upper()] for seq in seqs))]
        for i in range(len(seqs[0]))
    )


def bad_bases(seq):
    """Characters in seq that are neither a base nor an IUPAC ambiguity code."""
    return sorted({c for c in seq.upper() if c not in VALID_BASES})


def by_length(seqs):
    """Group sequences by length, in order of first appearance."""
    groups = {}
    for seq in seqs:
        groups.setdefault(len(seq), []).append(seq)
    return groups


def is_fasta(path):
    """True when the first non-blank line starts with '>'."""
    with open(path) as handle:
        for line in handle:
            if line.strip():
                return line.startswith(">")
    return False


def resolve_inputs(input_path):
    """Expand --input into the list of primer FASTAs to parse."""
    input_path = Path(input_path)
    if not input_path.is_dir():
        return [input_path]

    paths = sorted(p for p in input_path.iterdir() if p.suffix.lower() in FASTA_SUFFIXES)
    if not paths:
        sys.exit(
            f"No primer FASTA files ({', '.join(FASTA_SUFFIXES)}) found in {input_path} - "
            "--input must be a primer FASTA, a directory of primer FASTAs, or a samplesheet."
        )
    return paths


def group_by_prefix(records, path):
    """Sort validated records into {prefix: {'fwd': [...], 'rev': [...]}}, or report why not."""
    tagged = [(header, seq, split_header(header)) for header, seq in records]
    untagged = [header for header, _seq, (_prefix, direction) in tagged if direction is None]

    if untagged and not (len(untagged) == len(records) == 2):
        return None, (
            f"{path.name}: cannot tell which direction these records belong to - {untagged}. "
            "Add an fwd/forward/rev/reverse token to the header, or reduce the file to a "
            "single forward/reverse pair."
        )

    groups = {}
    if untagged:
        # No tokens anywhere and exactly one record per direction - unambiguous.
        groups[""] = {"fwd": [records[0][1]], "rev": [records[1][1]]}
        return groups, None

    for _header, seq, (prefix, direction) in tagged:
        groups.setdefault(prefix, {"fwd": [], "rev": []})[direction].append(seq)

    for prefix, directions in groups.items():
        if not directions["fwd"] or not directions["rev"]:
            return None, (
                f"{path.name}: prefix '{prefix}' has {len(directions['fwd'])} fwd and "
                f"{len(directions['rev'])} rev - a primer set needs both."
            )
    return groups, None


def rows_from_file(path, min_len, max_len):
    """Parse one primer FASTA into (rows, warnings, errors)."""
    records = list(read_fasta(path))
    if not records:
        return [], [], [f"{path.name}: no FASTA records found - is this really a FASTA file?"]

    for header, seq in records:
        if not seq:
            return [], [], [f"{path.name}: record '{header}' has an empty sequence."]
        bad = bad_bases(seq)
        if bad:
            return [], [], [
                f"{path.name}: record '{header}' contains non-nucleotide character(s) "
                f"{', '.join(bad)} - only ACGT/U and IUPAC ambiguity codes are allowed."
            ]

    groups, error = group_by_prefix(records, path)
    if error:
        return [], [], [error]

    rows = []
    warnings = []
    stem = UNSAFE_NAME_RE.sub("_", path.stem)
    multiple_prefixes = len(groups) > 1

    for prefix, directions in groups.items():
        fwd_groups = by_length(directions["fwd"])
        rev_groups = by_length(directions["rev"])
        combinations = [
            (collapse(fwd), collapse(rev))
            for fwd in fwd_groups.values()
            for rev in rev_groups.values()
        ]

        base = f"{stem}_{UNSAFE_NAME_RE.sub('_', prefix)}" if multiple_prefixes and prefix else stem
        if len(combinations) > 1:
            warnings.append(
                f"{path.name}: prefix '{prefix}' could not be collapsed - fwd lengths "
                f"{sorted(fwd_groups)}, rev lengths {sorted(rev_groups)}. Benchmarking "
                f"{len(combinations)} combinations separately."
            )

        for index, (fwd, rev) in enumerate(combinations, start=1):
            name = f"{base}_{index}" if len(combinations) > 1 else base
            rows.append({"primer": name, "fwd": fwd, "rev": rev, "min": min_len, "max": max_len})

    return rows, warnings, []


def collect_rows(paths, min_len, max_len):
    """Parse every primer FASTA, collecting all problems rather than stopping at the first."""
    rows, warnings, errors = [], [], []
    for path in paths:
        file_rows, file_warnings, file_errors = rows_from_file(Path(path), min_len, max_len)
        rows.extend(file_rows)
        warnings.extend(file_warnings)
        errors.extend(file_errors)

    # Nothing is benchmarked while any file is unusable, so a typo can never
    # silently drop a primer from the comparison.
    return ([], warnings, errors) if errors else (rows, warnings, errors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Primer FASTA, or directory of primer FASTAs")
    parser.add_argument("--min", required=True, type=int, help="Global minimum amplicon length")
    parser.add_argument("--max", required=True, type=int, help="Global maximum amplicon length")
    parser.add_argument("--out", required=True, help="Output samplesheet (TSV)")
    parser.add_argument("--warnings", help="Optional file to write parser warnings to")
    args = parser.parse_args()

    paths = resolve_inputs(args.input)
    rows, warnings, errors = collect_rows(paths, args.min, args.max)

    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if args.warnings:
        Path(args.warnings).write_text("".join(f"{warning}\n" for warning in warnings))

    if errors:
        detail = "".join(f"\n  {error}" for error in errors)
        sys.exit(
            f"{len(errors)} of {len(paths)} primer FASTA file(s) could not be used:{detail}\n\n"
            "No primer sets were run. Fix these and rerun."
        )

    with open(args.out, "w") as handle:
        handle.write("\t".join(TSV_HEADER) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in TSV_HEADER) + "\n")


if __name__ == "__main__":
    main()
