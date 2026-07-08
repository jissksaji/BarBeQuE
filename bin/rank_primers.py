#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import sys
import argparse

parser = argparse.ArgumentParser(description="Rank primers based on consensus files and teeliste targets.")
parser.add_argument(
    "-d", "--data-dir",
    type=str,
    default="/home/saj/jiss/barbeque/analysis",
    help="Path to the analysis data directory containing 'teeliste.tsv' and the 'consensus' folder."
)
args = parser.parse_args()

DATA_DIR = Path(args.data_dir)
teeliste_path = DATA_DIR / "teeliste.tsv"
consensus_dir = DATA_DIR / "consensus"

if not teeliste_path.exists() or not consensus_dir.exists():
    print(f"Could not find the analysis data directories in {DATA_DIR}.")
    sys.exit(1)

print("Loading teeliste targets...")
teeliste_df = pd.read_csv(teeliste_path, sep="\t", header=0)
teeliste_df.columns = ["german_name", "latin_name", "taxid"]
teeliste_taxids = set(teeliste_df["taxid"].astype(str).str.strip())
TOTAL_TARGETS = len(teeliste_taxids)

print(f"Loaded {TOTAL_TARGETS} targets from teeliste.tsv.")

COLUMNS = [
    "cluster_id", "accession", "accession_taxid", "accession_name",
    "assigned_name", "assigned_taxid", "assigned_rank", "disambiguation"
]

results = []

consensus_files = list(consensus_dir.glob("*_consensus.tsv"))
print(f"Found {len(consensus_files)} consensus files. Analyzing...")

for file_path in consensus_files:
    primer_name = file_path.name.replace("_consensus.tsv", "")
    
    try:
        df = pd.read_csv(file_path, sep="\t", header=None, names=COLUMNS, on_bad_lines="skip")
    except Exception as e:
        print(f"Error reading {file_path.name}: {e}")
        continue
        
    if df.empty:
        continue

    # 1. Calculate Resolution Efficiency (% of clusters resolving to species)
    # Group by cluster_id to get one row per cluster
    clusters_df = df[["cluster_id", "assigned_rank"]].drop_duplicates(subset=["cluster_id"])
    total_clusters = len(clusters_df)
    
    if total_clusters == 0:
        continue
        
    species_clusters = len(clusters_df[clusters_df["assigned_rank"].str.strip() == "species"])
    resolution_pct = (species_clusters / total_clusters) * 100
    
    genus_clusters = len(clusters_df[clusters_df["assigned_rank"].str.strip() == "genus"])
    resolution_genus_pct = (genus_clusters / total_clusters) * 100

    family_clusters = len(clusters_df[clusters_df["assigned_rank"].str.strip() == "family"])
    resolution_family_pct = (family_clusters / total_clusters) * 100

    # 2. Calculate Coverage (% of teeliste targets amplified)
    # A target is "amplified" if its taxid is found in accession_taxid OR assigned_taxid
    accession_taxids = set(df["accession_taxid"].astype(str).str.strip().unique())
    assigned_taxids = set(df["assigned_taxid"].astype(str).str.strip().unique())
    
    all_amplified_taxids = accession_taxids.union(assigned_taxids)
    total_amplified_taxids_count = len(all_amplified_taxids)
    
    amplified_count = 0
    for tid in teeliste_taxids:
        if tid in accession_taxids or tid in assigned_taxids:
            amplified_count += 1
            
    coverage_pct = (amplified_count / TOTAL_TARGETS) * 100

    # 3. Calculate Teeliste-Specific Resolution
    teeliste_subset = df[df["accession_taxid"].astype(str).str.strip().isin(teeliste_taxids) | df["assigned_taxid"].astype(str).str.strip().isin(teeliste_taxids)]
    teeliste_clusters_df = teeliste_subset[["cluster_id", "assigned_rank"]].drop_duplicates(subset=["cluster_id"])
    total_teeliste_clusters = len(teeliste_clusters_df)
    
    if total_teeliste_clusters > 0:
        teeliste_species = len(teeliste_clusters_df[teeliste_clusters_df["assigned_rank"].str.strip() == "species"])
        teeliste_genus = len(teeliste_clusters_df[teeliste_clusters_df["assigned_rank"].str.strip() == "genus"])
        teeliste_family = len(teeliste_clusters_df[teeliste_clusters_df["assigned_rank"].str.strip() == "family"])
        teeliste_res_species_pct = (teeliste_species / total_teeliste_clusters) * 100
        teeliste_res_genus_pct = (teeliste_genus / total_teeliste_clusters) * 100
        teeliste_res_family_pct = (teeliste_family / total_teeliste_clusters) * 100
    else:
        teeliste_res_species_pct = 0.0
        teeliste_res_genus_pct = 0.0
        teeliste_res_family_pct = 0.0

    # 4. Calculate Combined Score (Simple Average of Species Res and Target Cov)
    combined_score = (resolution_pct + coverage_pct) / 2

    results.append({
        "Primer Name": primer_name,
        "Total Clusters": total_clusters,
        "Amplified TaxIDs": total_amplified_taxids_count,
        "Resolution (Species %)": round(resolution_pct, 2),
        "Resolution (Genus %)": round(resolution_genus_pct, 2),
        "Resolution (Family %)": round(resolution_family_pct, 2),
        "Coverage (Teeliste %)": round(coverage_pct, 2),
        "Teeliste Res (Species %)": round(teeliste_res_species_pct, 2),
        "Teeliste Res (Genus %)": round(teeliste_res_genus_pct, 2),
        "Teeliste Res (Family %)": round(teeliste_res_family_pct, 2),
        "Combined Score": round(combined_score, 2)
    })

if not results:
    print("No valid primer data could be parsed.")
    sys.exit(1)

# Create DataFrame and sort by Combined Score descending
results_df = pd.DataFrame(results).sort_values(by="Combined Score", ascending=False).reset_index(drop=True)

# Format the percentages for pretty printing
results_df["Resolution (Species %)"] = results_df["Resolution (Species %)"].astype(str) + "%"
results_df["Resolution (Genus %)"] = results_df["Resolution (Genus %)"].astype(str) + "%"
results_df["Resolution (Family %)"] = results_df["Resolution (Family %)"].astype(str) + "%"
results_df["Coverage (Teeliste %)"] = results_df["Coverage (Teeliste %)"].astype(str) + "%"
results_df["Teeliste Res (Species %)"] = results_df["Teeliste Res (Species %)"].astype(str) + "%"
results_df["Teeliste Res (Genus %)"] = results_df["Teeliste Res (Genus %)"].astype(str) + "%"
results_df["Teeliste Res (Family %)"] = results_df["Teeliste Res (Family %)"].astype(str) + "%"

print("\n=== TOP 20 PRIMERS BY COMBINED SCORE ===")
print(results_df.head(20).to_string(index=False))

# Save to file
output_file = DATA_DIR / "primer_rankings.tsv"
results_df.to_csv(output_file, sep="\t", index=False)
print(f"\nSaved full rankings for all {len(results_df)} primers to: {output_file}")
