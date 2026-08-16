#!/usr/bin/env python3
"""Remove blocklisted taxids from a FASTA database."""

import argparse
import re


def strip_version(accession):
    """Turn 'MK123456.1' or a SINTAX id into the base accession."""
    accession = accession.split(";", 1)[0]
    return re.sub(r"(\.[0-9]+)+$", "", accession)


def load_blocked_taxids(path):
    """Read one taxid per line. Blank lines and # comments are ignored."""
    blocked_taxids = set()
    with open(path) as handle:
        for line in handle:
            taxid = line.split("#", 1)[0].strip()
            if taxid:
                blocked_taxids.add(taxid)
    return blocked_taxids


def load_fasta_accessions(path):
    """Read the accessions that are actually present in the selected FASTA."""
    accessions = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                accessions.add(strip_version(line[1:].split()[0]))
    return accessions


def load_blocked_accessions(path, blocked_taxids, database_accessions=None):
    """Find accessions assigned to a blocked taxid.

    Both two-column files (accession, taxid) and standard four-column NCBI
    accession2taxid files are supported.
    """
    blocked_accessions = set()
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 3:
                accession, taxid = fields[0], fields[2]
            elif len(fields) >= 2:
                accession, taxid = fields[0], fields[1]
            else:
                continue

            accession = strip_version(accession)
            is_in_database = (
                database_accessions is None or accession in database_accessions
            )
            if taxid in blocked_taxids and is_in_database:
                blocked_accessions.add(accession)

    return blocked_accessions


def filter_fasta(input_path, output_path, blocked_accessions):
    """Copy unblocked FASTA records to output_path and return record counts."""
    total = 0
    removed = 0
    record = []
    is_blocked = False

    with open(input_path) as source, open(output_path, "w") as output:
        for line in source:
            if line.startswith(">"):
                if record and not is_blocked:
                    output.writelines(record)

                total += 1
                accession = strip_version(line[1:].split()[0])
                is_blocked = accession in blocked_accessions
                if is_blocked:
                    removed += 1
                record = [line]
            else:
                record.append(line)

        if record and not is_blocked:
            output.writelines(record)

    return total, removed


def write_summary(path, total, removed, blocked_taxids, blocked_accessions):
    """Write a small, human-readable TSV summary."""
    with open(path, "w") as output:
        output.write("metric\tcount\n")
        output.write(f"input_records\t{total}\n")
        output.write(f"kept_records\t{total - removed}\n")
        output.write(f"removed_records\t{removed}\n")
        output.write(f"blocklist_taxids\t{len(blocked_taxids)}\n")
        output.write(f"blocked_accessions\t{len(blocked_accessions)}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="Input FASTA database")
    parser.add_argument(
        "--accession-taxid",
        required=True,
        help="Two- or four-column accession-to-taxid file",
    )
    parser.add_argument("--blocklist", required=True, help="FooDMe2 taxid blocklist")
    parser.add_argument("--output", required=True, help="Filtered FASTA output")
    parser.add_argument("--summary", required=True, help="TSV summary output")
    args = parser.parse_args()

    blocked_taxids = load_blocked_taxids(args.blocklist)
    database_accessions = load_fasta_accessions(args.fasta)
    blocked_accessions = load_blocked_accessions(
        args.accession_taxid, blocked_taxids, database_accessions
    )
    total, removed = filter_fasta(args.fasta, args.output, blocked_accessions)
    write_summary(
        args.summary, total, removed, blocked_taxids, blocked_accessions
    )


if __name__ == "__main__":
    main()
