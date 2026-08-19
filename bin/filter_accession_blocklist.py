#!/usr/bin/env python3
"""Remove listed accessions from parsed OBI-PCR results and amplicon FASTA.

The accession blocklist contains one accession per line. Blank lines and text
after ``#`` are ignored. Matching is case- and version-insensitive, so an entry
such as ``MK123456`` also removes ``MK123456.1``.
"""

import argparse
import re


def normalize_accession(accession):
    """Return a case- and version-normalized accession identifier."""
    # OBI-PCR FASTA headers may contain a description after whitespace or
    # semicolon-delimited annotations. Only the leading identifier is relevant.
    accession = accession.strip().split()[0].split(";", 1)[0]

    # Some inputs contain compound numeric suffixes (for example
    # AY846379.1.1791). Removing every trailing numeric component makes a
    # versionless blocklist entry match all versions of the same accession.
    accession = re.sub(r"(\.[0-9]+)+$", "", accession)
    return accession.upper()


def load_accession_blocklist(path):
    """Load unique accessions, ignoring blank lines and comments."""
    accessions = set()
    with open(path) as handle:
        for line in handle:
            token = line.split("#", 1)[0].strip()
            if token:
                accessions.add(normalize_accession(token))
    return accessions


def filter_fasta(input_path, output_path, blocked_accessions):
    """Filter FASTA records and return counts plus matched accessions."""
    total = 0
    removed = 0
    matched = set()
    record = []
    is_blocked = False

    # Stream records instead of parsing the whole FASTA into memory. Keeping the
    # original lines also preserves wrapping and header descriptions verbatim.
    with open(input_path) as source, open(output_path, "w") as output:
        for line in source:
            if line.startswith(">"):
                # A new header closes the preceding record, so write that record
                # only after its blocked status is known.
                if record and not is_blocked:
                    output.writelines(record)

                total += 1
                accession = normalize_accession(line[1:])
                is_blocked = accession in blocked_accessions
                if is_blocked:
                    removed += 1
                    matched.add(accession)
                record = [line]
            else:
                record.append(line)

        if record and not is_blocked:
            output.writelines(record)

    return total, removed, matched


def filter_tsv(input_path, output_path, blocked_accessions):
    """Filter parsed OBI-PCR rows by the Sequence_ID first column."""
    total = 0
    removed = 0
    matched = set()

    with open(input_path) as source, open(output_path, "w") as output:
        header = source.readline()
        if not header:
            raise ValueError(f"Parsed OBI-PCR table is empty: {input_path}")
        if header.rstrip("\r\n").split("\t", 1)[0] != "Sequence_ID":
            raise ValueError(
                f"Parsed OBI-PCR table must start with a Sequence_ID column: {input_path}"
            )
        output.write(header)

        # Inspect only Sequence_ID and otherwise copy the source row unchanged;
        # using a table library here would add a dependency without simplifying
        # this fixed first-column filter.
        for line in source:
            if not line.strip():
                continue
            total += 1
            accession = normalize_accession(line.split("\t", 1)[0])
            if accession in blocked_accessions:
                removed += 1
                matched.add(accession)
            else:
                output.write(line)

    return total, removed, matched


def write_summary(
    path,
    listed_accessions,
    fasta_counts,
    tsv_counts,
    matched_accessions,
):
    """Write removal counts for auditing the post-OBI-PCR filter."""
    fasta_total, fasta_removed = fasta_counts
    tsv_total, tsv_removed = tsv_counts
    with open(path, "w") as output:
        output.write("metric\tcount\n")
        output.write(f"listed_accessions\t{len(listed_accessions)}\n")
        output.write(f"matched_accessions\t{len(matched_accessions)}\n")
        output.write(
            f"unmatched_accessions\t{len(listed_accessions - matched_accessions)}\n"
        )
        output.write(f"fasta_input_records\t{fasta_total}\n")
        output.write(f"fasta_kept_records\t{fasta_total - fasta_removed}\n")
        output.write(f"fasta_removed_records\t{fasta_removed}\n")
        output.write(f"parsed_input_rows\t{tsv_total}\n")
        output.write(f"parsed_kept_rows\t{tsv_total - tsv_removed}\n")
        output.write(f"parsed_removed_rows\t{tsv_removed}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="OBI-PCR amplicon FASTA")
    parser.add_argument("--tsv", required=True, help="Parsed OBI-PCR TSV")
    parser.add_argument(
        "--accession-blocklist",
        required=True,
        help="One accession per line; blank lines and # comments are ignored",
    )
    parser.add_argument("--fasta-output", required=True)
    parser.add_argument("--tsv-output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    blocked_accessions = load_accession_blocklist(args.accession_blocklist)
    if not blocked_accessions:
        parser.error("the accession blocklist contains no accessions")

    fasta_total, fasta_removed, fasta_matched = filter_fasta(
        args.fasta,
        args.fasta_output,
        blocked_accessions,
    )
    tsv_total, tsv_removed, tsv_matched = filter_tsv(
        args.tsv,
        args.tsv_output,
        blocked_accessions,
    )
    write_summary(
        args.summary,
        blocked_accessions,
        (fasta_total, fasta_removed),
        (tsv_total, tsv_removed),
        # An accession counts as matched if it appeared in either representation.
        # Separate FASTA/TSV removal counts still expose any disagreement.
        fasta_matched | tsv_matched,
    )


if __name__ == "__main__":
    main()
