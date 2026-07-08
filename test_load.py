import pandas as pd
from pathlib import Path
import sys

DATA_DIR = Path("/home/saj/jiss/barbeque/results_teupdated")

COLUMNS = [
    "cluster_id",
    "accession",
    "accession_taxid",
    "accession_name",
    "assigned_name",
    "assigned_taxid",
    "assigned_rank",
    "disambiguation"
]

DTYPES = {
    "cluster_id": "Int64",
    "accession": "string",
    "accession_taxid": "Int64",
    "accession_name": "string",
    "assigned_name": "string",
    "assigned_taxid": "string",
    "assigned_rank": "string",
    "disambiguation": "string"
}

consensus_dir = DATA_DIR / "consensus"
print("Reading from:", consensus_dir)

dfs = []
for file in sorted(consensus_dir.glob("*.cluster_consensus.tsv")):
    print("Reading", file)
    tmp = pd.read_csv(
        file,
        sep="\t",
        header=None,
        names=COLUMNS,
        dtype=DTYPES
    )
    tmp["accession_name"] = "Unknown"
    tmp["primer"] = Path(file.stem).stem
    dfs.append(tmp)

df = pd.concat(dfs, ignore_index=True)
print("Dataframe shape:", df.shape)

length_dir = DATA_DIR / "amplicon_lengths"
dfs_len = []
for file in sorted(length_dir.glob("*.tsv")):
    tmp = pd.read_csv(file, sep="\t", header=None, names=["header", "length"])
    stem = file.stem.replace(".amplicon_lengths", "")
    tmp["primer"] = stem
    dfs_len.append(tmp)

filtered_len = pd.concat(dfs_len, ignore_index=True)

merged_len = pd.merge(filtered_len, df[["accession", "assigned_name", "assigned_rank"]].drop_duplicates(), left_on="header", right_on="accession", how="inner")
print("Merged shape:", merged_len.shape)
