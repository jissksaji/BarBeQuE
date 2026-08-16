#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description="Plot amplicons and multi-amplicons directly from parsed obipcr output.")
    parser.add_argument("-d", "--data-dir", default="BarBeQuE/results/parsed_obipcr", help="Path to parsed_obipcr directory")
    parser.add_argument("-o", "--output", default="BarBeQuE/results/obipcr_amplicons_plot.png", help="Output plot filename")
    parser.add_argument("-n", "--top-n", type=int, default=30, help="Number of top primers to plot")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        print(f"Error: Directory {data_dir} not found.")
        return

    results = []
    print("Scanning parsed_obipcr files...")
    for file in data_dir.glob("*.tsv"):
        try:
            # Read only Sequence_ID to save memory
            df = pd.read_csv(file, sep="\t", usecols=["Sequence_ID"])
            total_amplicons = len(df)
            unique_accessions = df["Sequence_ID"].nunique()
            multi_amplicons = total_amplicons - unique_accessions
            
            # Primer name (strip _euphyllophyta if present for cleaner labels)
            primer_name = file.stem.replace("_euphyllophyta", "")
            
            results.append({
                "Primer": primer_name,
                "Unique Accessions": unique_accessions,
                "Multi-Amplicons": multi_amplicons,
                "Total Amplicons": total_amplicons
            })
        except Exception as e:
            print(f"Skipping {file.name}: {e}")

    if not results:
        print("No valid data found to plot.")
        return

    df_results = pd.DataFrame(results)
    # Sort by total amplicons descending
    df_results = df_results.sort_values(by="Total Amplicons", ascending=False).head(args.top_n)

    # Plotting
    plt.figure(figsize=(14, 8))
    
    # Create simple bar chart
    bars = plt.bar(df_results["Primer"], df_results["Total Amplicons"], color='#1f77b4', label='Total Amplicons')
    
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Total Number of Amplicons')
    plt.title(f'Top {args.top_n} Primers by Total Amplicons (Raw obipcr)')
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(args.output, dpi=300)
    print(f"Plot successfully saved to {args.output}")

if __name__ == "__main__":
    main()
