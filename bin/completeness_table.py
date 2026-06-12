#!/usr/bin/env python3
"""
completeness_table.py

Takes a list of species of interest and the taxids_counts file from BUILD_DB_TAXIDS.
Outputs a 10-column TSV describing how well each species' genus and family
are represented in the reference database.

Usage:
    python completeness_table.py \
        --species "Camellia sinensis+Camellia japonica" \
        --taxids_counts db_taxids_counts.tsv \
        --taxdump /path/to/new_taxdump \
        --output completeness.tsv

    or with a .txt file:
    python completeness_table.py \
        --species species_list.txt \
        --taxids_counts db_taxids_counts.tsv \
        --taxdump /path/to/new_taxdump \
        --output completeness.tsv
"""

import argparse
import os
import tarfile
import tempfile
from ete3 import NCBITaxa


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a completeness table for species of interest")

    parser.add_argument("--species", required=True,
        help="Species names separated by '+', or path to a .txt file with one species per line")

    parser.add_argument("--taxids_counts", required=True,
        help="TSV file from BUILD_DB_TAXIDS with columns: taxid, count")

    parser.add_argument("--taxdump", required=True,
        help="Path to directory containing names.dmp and nodes.dmp")

    parser.add_argument("--output", required=True,
        help="Output TSV file path")

    return parser.parse_args()


def load_species_list(species_input):
    # if it looks like a file path and the file exists, read one species per line
    if os.path.isfile(species_input):
        with open(species_input) as file:
            species_list = [line.strip() for line in file if line.strip()]
    else:
        # otherwise split on + separator
        species_list = [s.strip() for s in species_input.split("+") if s.strip()]

    return species_list


def load_taxids_counts(taxids_counts_path):
    # read the taxids_counts file into a dict: taxid (int) -> sequence count (int)
    taxid_to_count = {}

    with open(taxids_counts_path) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            taxid = int(parts[0])
            count = int(parts[1])
            taxid_to_count[taxid] = count

    return taxid_to_count


def get_rank_taxid(taxid, target_rank, ncbi):
    # walk up the lineage of a taxid and return the taxid at the target rank
    # returns None if the rank is not found
    lineage = ncbi.get_lineage(taxid)
    rank_map = ncbi.get_rank(lineage)

    for ancestor_taxid in lineage:
        if rank_map.get(ancestor_taxid) == target_rank:
            return ancestor_taxid

    return None


def get_species_under_taxid(parent_taxid, ncbi):
    # get all taxids that are descendants of parent_taxid and have rank 'species'
    all_descendants = ncbi.get_descendant_taxa(parent_taxid, collapse_subspecies=False)
    rank_map = ncbi.get_rank(all_descendants)

    species_taxids = [t for t in all_descendants if rank_map.get(t) == "species"]
    return species_taxids


def get_name(taxid, ncbi):
    # return the scientific name for a taxid
    name_map = ncbi.get_taxid_translator([taxid])
    return name_map.get(taxid, str(taxid))


def build_row(species_name, taxid_to_count, ncbi):
    # resolve the species name to a taxid
    name_to_taxid = ncbi.get_name_translator([species_name])

    if species_name not in name_to_taxid:
        print(f"  WARNING: could not resolve '{species_name}' in NCBI taxonomy, skipping")
        return None

    species_taxid = name_to_taxid[species_name][0]

    # get the genus and family taxids by walking up the lineage
    genus_taxid = get_rank_taxid(species_taxid, "genus", ncbi)
    family_taxid = get_rank_taxid(species_taxid, "family", ncbi)


    # column 2: number of sequences for this species in the DB
    seq_count = taxid_to_count.get(species_taxid, 0)


    # get all species taxids in the DB (any taxid that has a count)
    db_taxids = set(taxid_to_count.keys())


    # columns 3, 9: species in DB sharing the same genus
    if genus_taxid is not None:
        all_genus_species = get_species_under_taxid(genus_taxid, ncbi)
        db_genus_species = [t for t in all_genus_species if t in db_taxids]
        db_genus_count = len(db_genus_species)
        ncbi_genus_count = len(all_genus_species)
        db_genus_names = "; ".join(sorted(get_name(t, ncbi) for t in db_genus_species))
    else:
        db_genus_count = 0
        ncbi_genus_count = 0
        db_genus_names = ""

    # column 5: percentage of genus in DB
    if ncbi_genus_count > 0:
        genus_percent = round(db_genus_count / ncbi_genus_count * 100, 2)
    else:
        genus_percent = 0.0


    # columns 6, 10: species in DB sharing the same family
    if family_taxid is not None:
        all_family_species = get_species_under_taxid(family_taxid, ncbi)
        db_family_species = [t for t in all_family_species if t in db_taxids]
        db_family_count = len(db_family_species)
        ncbi_family_count = len(all_family_species)
        db_family_names = "; ".join(sorted(get_name(t, ncbi) for t in db_family_species))
    else:
        db_family_count = 0
        ncbi_family_count = 0
        db_family_names = ""

    # column 8: percentage of family in DB
    if ncbi_family_count > 0:
        family_percent = round(db_family_count / ncbi_family_count * 100, 2)
    else:
        family_percent = 0.0


    # assemble the row in column order 1-10
    row = [
        species_name,           # col 1: species name
        seq_count,              # col 2: sequences in DB for this species
        db_genus_count,         # col 3: species in DB sharing genus
        ncbi_genus_count,       # col 4: total species in genus (NCBI)
        genus_percent,          # col 5: % of genus in DB
        db_family_count,        # col 6: species in DB sharing family
        ncbi_family_count,      # col 7: total species in family (NCBI)
        family_percent,         # col 8: % of family in DB
        db_genus_names,         # col 9: names of genus species in DB
        db_family_names,        # col 10: names of family species in DB
    ]

    return row


def main():
    args = parse_args()

    print("Loading NCBI taxonomy...")
    
    taxdump_dir = args.taxdump
    db_path = "taxa.sqlite"
    tar_gz_path = os.path.join(taxdump_dir, "taxdump.tar.gz")
    tmp_tar_path = None
    
    if not os.path.isfile(tar_gz_path):
        names_path = os.path.join(taxdump_dir, "names.dmp")
        if os.path.isfile(names_path):
            print("  Packaging .dmp files into a temporary tar.gz...")
            fd, tmp_tar_path = tempfile.mkstemp(suffix=".tar.gz", dir=".")
            os.close(fd)
            with tarfile.open(tmp_tar_path, "w:gz") as tar:
                for fname in os.listdir(taxdump_dir):
                    if fname.endswith(".dmp"):
                        tar.add(os.path.join(taxdump_dir, fname), arcname=fname)
            tar_gz_path = tmp_tar_path
        else:
            tar_gz_path = None

    print("  Initializing NCBITaxa...")
    if tar_gz_path:
        ncbi = NCBITaxa(dbfile=db_path, taxdump_file=tar_gz_path)
    else:
        ncbi = NCBITaxa()
        
    if tmp_tar_path and os.path.exists(tmp_tar_path):
        os.remove(tmp_tar_path)

    print(f"Loading taxids_counts from {args.taxids_counts}...")
    taxid_to_count = load_taxids_counts(args.taxids_counts)
    print(f"  {len(taxid_to_count)} taxids loaded")

    print(f"Loading species list...")
    species_list = load_species_list(args.species)
    print(f"  {len(species_list)} species to process")

    header = [
        "species_name",
        "db_seq_count",
        "db_genus_species_count",
        "ncbi_genus_species_count",
        "genus_coverage_pct",
        "db_family_species_count",
        "ncbi_family_species_count",
        "family_coverage_pct",
        "db_genus_species_names",
        "db_family_species_names",
    ]

    rows = []
    for species_name in species_list:
        print(f"  Processing: {species_name}")
        row = build_row(species_name, taxid_to_count, ncbi)
        if row is not None:
            rows.append(row)

    print(f"Writing output to {args.output}...")
    with open(args.output, "w") as out_file:
        out_file.write("\t".join(header) + "\n")
        for row in rows:
            out_file.write("\t".join(str(x) for x in row) + "\n")

    print(f"Done. {len(rows)} rows written.")


if __name__ == "__main__":
    main()