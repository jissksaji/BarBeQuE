#!/usr/bin/env python3
"""Plot Teeliste genus-level coverage against database representation.

The output is a standalone interactive HTML figure; Streamlit is not needed.
Coverage is matched by exact TaxID using the BarBeQuE consensus output.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


CONSENSUS_COLUMNS = [
    "cluster_id",
    "accession",
    "accession_taxid",
    "accession_name",
    "assigned_name",
    "assigned_taxid",
    "assigned_rank",
    "disambiguation",
]

GENUS_OR_FINER = {
    "genus",
    "subgenus",
    "species",
    "subspecies",
    "variety",
    "subvariety",
    "forma",
    "form",
    "strain",
}

STATUS_ORDER = [
    "Amplified — genus or finer",
    "Amplified — broader/unresolved",
    "In database, not amplified",
    "Missing from database",
]

STATUS_COLORS = {
    "Amplified — genus or finer": "#1b7837",
    "Amplified — broader/unresolved": "#80cdc1",
    "In database, not amplified": "#ff7f0e",
    "Missing from database": "#bdbdbd",
}


def clean_taxid(series):
    return series.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)


def is_genus_or_finer(ranks):
    observed = {
        str(rank).strip().lower()
        for rank in ranks
        if pd.notna(rank) and str(rank).strip()
    }
    return bool(observed & GENUS_OR_FINER)


def build_plot_data(teeliste_path, db_counts_path, consensus_path):
    targets = pd.read_csv(teeliste_path, sep="\t", header=0).iloc[:, :3].copy()
    targets.columns = ["german_name", "latin_name", "taxid"]
    targets["taxid"] = clean_taxid(targets["taxid"])
    targets = targets.drop_duplicates("taxid").reset_index(drop=True)

    db_counts = pd.read_csv(
        db_counts_path,
        sep="\t",
        header=None,
        names=["taxid", "database_sequences"],
        usecols=[0, 1],
    )
    db_counts["taxid"] = clean_taxid(db_counts["taxid"])
    db_counts["database_sequences"] = pd.to_numeric(
        db_counts["database_sequences"], errors="coerce"
    ).fillna(0)
    db_count_map = db_counts.groupby("taxid")["database_sequences"].sum()

    consensus = pd.read_csv(
        consensus_path,
        sep="\t",
        header=None,
        names=CONSENSUS_COLUMNS,
        dtype="string",
        on_bad_lines="skip",
    )
    consensus["accession_taxid"] = clean_taxid(consensus["accession_taxid"])
    consensus["assigned_taxid"] = clean_taxid(consensus["assigned_taxid"])

    rows = []
    for target in targets.itertuples(index=False):
        direct_matches = consensus[consensus["accession_taxid"].eq(target.taxid)]
        coverage_matches = consensus[
            consensus["accession_taxid"].eq(target.taxid)
            | consensus["assigned_taxid"].eq(target.taxid)
        ]
        # Keep the two sequence bars comparable: both count sequences assigned
        # directly to the target TaxID. Assigned-LCA matches still contribute
        # to target coverage and resolution, but not to this sequence count.
        amplified_sequences = direct_matches["accession"].dropna().nunique()
        amplified = not coverage_matches.empty
        genus_or_finer = amplified and is_genus_or_finer(
            coverage_matches["assigned_rank"]
        )
        database_sequences = int(db_count_map.get(target.taxid, 0))

        if genus_or_finer:
            status = STATUS_ORDER[0]
        elif amplified:
            status = STATUS_ORDER[1]
        elif database_sequences > 0:
            status = STATUS_ORDER[2]
        else:
            status = STATUS_ORDER[3]

        ranks = sorted(
            {
                str(rank).strip()
                for rank in coverage_matches["assigned_rank"]
                if pd.notna(rank) and str(rank).strip()
            }
        )
        rows.append(
            {
                "taxid": target.taxid,
                "german_name": target.german_name,
                "latin_name": target.latin_name,
                "database_sequences": database_sequences,
                "amplified_sequences": amplified_sequences,
                "best_observed_ranks": "; ".join(ranks) or "not resolved",
                "status": status,
            }
        )

    result = pd.DataFrame(rows)
    result["database_log"] = np.log10(result["database_sequences"] + 1)
    result["amplified_log"] = np.log10(result["amplified_sequences"] + 1)
    return result


def make_figure(plot_data, title):
    total = len(plot_data)
    genus_count = int((plot_data["status"] == STATUS_ORDER[0]).sum())
    coverage = genus_count / total * 100 if total else 0.0

    figure = px.scatter(
        plot_data,
        x="database_log",
        y="amplified_log",
        color="status",
        color_discrete_map=STATUS_COLORS,
        category_orders={"status": STATUS_ORDER},
        hover_name="latin_name",
        hover_data={
            "german_name": True,
            "taxid": True,
            "database_sequences": ":,",
            "amplified_sequences": ":,",
            "best_observed_ranks": True,
            "database_log": False,
            "amplified_log": False,
        },
        labels={
            "database_log": "Sequences in database (log10(count + 1))",
            "amplified_log": "Sequences amplified (log10(count + 1))",
            "status": "Target result",
            "best_observed_ranks": "Observed consensus ranks",
        },
        title=f"{title}<br><sup>Genus-or-finer coverage: {genus_count}/{total} ({coverage:.1f}%)</sup>",
    )
    figure.update_traces(marker={"size": 11, "opacity": 0.82})
    figure.update_layout(
        template="plotly_white",
        height=650,
        legend_title_text="Target result",
        margin={"l": 70, "r": 30, "t": 100, "b": 70},
    )
    return figure


def summarize_primers(teeliste_path, db_counts_path, consensus_dir, db_name=None):
    rows = []
    suffix = ".cluster_consensus.tsv"
    for consensus_path in sorted(consensus_dir.glob(f"*{suffix}")):
        primer = consensus_path.name[: -len(suffix)]
        if db_name and primer.endswith(f"_{db_name}"):
            primer = primer[: -(len(db_name) + 1)]

        data = build_plot_data(teeliste_path, db_counts_path, consensus_path)
        total_targets = len(data)
        covered_targets = int(
            data["status"].isin(STATUS_ORDER[:2]).sum()
        )
        genus_targets = int((data["status"] == STATUS_ORDER[0]).sum())
        rows.append(
            {
                "primer": primer,
                "database_sequences": int(data["database_sequences"].sum()),
                "amplified_sequences": int(data["amplified_sequences"].sum()),
                "covered_targets": covered_targets,
                "genus_or_finer_targets": genus_targets,
                "total_targets": total_targets,
                "genus_coverage_pct": (
                    genus_targets / total_targets * 100 if total_targets else 0.0
                ),
            }
        )

    if not rows:
        raise ValueError(f"No *{suffix} files found in {consensus_dir}")
    return pd.DataFrame(rows).sort_values(
        ["genus_coverage_pct", "amplified_sequences"],
        ascending=[True, True],
    )


def make_primer_comparison(summary, title):
    figure = make_subplots(
        rows=1,
        cols=2,
        shared_yaxes=True,
        horizontal_spacing=0.08,
        subplot_titles=(
            "Target-list sequence counts",
            "Targets resolved to genus or finer",
        ),
    )
    figure.add_trace(
        go.Bar(
            x=summary["database_sequences"],
            y=summary["primer"],
            name="Sequences in database",
            orientation="h",
            marker_color="#c7c7c7",
            hovertemplate="%{y}<br>Database sequences: %{x:,}<extra></extra>",
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Bar(
            x=summary["amplified_sequences"],
            y=summary["primer"],
            name="Sequences amplified",
            orientation="h",
            marker_color="#3182bd",
            hovertemplate="%{y}<br>Amplified sequences: %{x:,}<extra></extra>",
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Bar(
            x=summary["genus_coverage_pct"],
            y=summary["primer"],
            name="Genus-or-finer coverage",
            orientation="h",
            marker_color="#1b7837",
            text=[
                f"{covered}/{total} ({pct:.1f}%)"
                for covered, total, pct in zip(
                    summary["genus_or_finer_targets"],
                    summary["total_targets"],
                    summary["genus_coverage_pct"],
                )
            ],
            textposition="outside",
            customdata=np.column_stack(
                [summary["covered_targets"], summary["total_targets"]]
            ),
            hovertemplate=(
                "%{y}<br>Genus-or-finer: %{text}"
                "<br>Amplified targets at any rank: %{customdata[0]}/%{customdata[1]}"
                "<extra></extra>"
            ),
        ),
        row=1,
        col=2,
    )
    figure.update_layout(
        title=title,
        template="plotly_white",
        barmode="overlay",
        height=max(550, 48 * len(summary) + 180),
        margin={"l": 220, "r": 80, "t": 100, "b": 70},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.08},
    )
    figure.update_xaxes(title_text="Number of sequences", row=1, col=1)
    figure.update_xaxes(
        title_text="Teeliste coverage (%)", range=[0, 105], row=1, col=2
    )
    figure.update_yaxes(title_text="Primer", row=1, col=1)
    return figure


def write_png(plot_data, title, output_path):
    """Write a static version of the coverage scatter without Kaleido."""
    import matplotlib.pyplot as plt

    total = len(plot_data)
    genus_count = int((plot_data["status"] == STATUS_ORDER[0]).sum())
    coverage = genus_count / total * 100 if total else 0.0

    figure, axis = plt.subplots(figsize=(12, 8))
    for status in STATUS_ORDER:
        subset = plot_data[plot_data["status"] == status]
        if subset.empty:
            continue
        axis.scatter(
            subset["database_log"],
            subset["amplified_log"],
            label=f"{status} (n={len(subset)})",
            color=STATUS_COLORS[status],
            s=72,
            alpha=0.82,
            edgecolors="white",
            linewidths=0.6,
        )

    axis.set_title(
        f"{title}\nGenus-or-finer coverage: "
        f"{genus_count}/{total} ({coverage:.1f}%)",
        pad=16,
    )
    axis.set_xlabel("Sequences in database (log10(count + 1))")
    axis.set_ylabel("Sequences amplified (log10(count + 1))")
    axis.grid(alpha=0.2)
    axis.legend(title="Target result", loc="upper left", frameon=True)
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def write_primer_comparison_png(summary, title, output_path):
    """Write the all-primer comparison as a static two-panel PNG."""
    import matplotlib.pyplot as plt

    height = max(6.5, 0.52 * len(summary) + 2.2)
    figure, (count_axis, coverage_axis) = plt.subplots(
        1,
        2,
        figsize=(15, height),
        sharey=True,
        gridspec_kw={"width_ratios": [1.25, 1]},
    )
    positions = np.arange(len(summary))

    count_axis.barh(
        positions,
        summary["database_sequences"],
        color="#c7c7c7",
        label="Sequences in database",
    )
    count_axis.barh(
        positions,
        summary["amplified_sequences"],
        color="#3182bd",
        alpha=0.9,
        label="Sequences amplified",
    )
    count_axis.set_yticks(positions, labels=summary["primer"])
    count_axis.set_xlabel("Number of target-list sequences")
    count_axis.set_title("Database vs amplified sequences")
    count_axis.legend(loc="lower right")
    count_axis.grid(axis="x", alpha=0.2)

    coverage_axis.barh(
        positions,
        summary["genus_coverage_pct"],
        color="#1b7837",
    )
    for position, row in enumerate(summary.itertuples(index=False)):
        coverage_axis.text(
            min(row.genus_coverage_pct + 1.2, 96),
            position,
            f"{row.genus_or_finer_targets}/{row.total_targets} "
            f"({row.genus_coverage_pct:.1f}%)",
            va="center",
            fontsize=9,
        )
    coverage_axis.set_xlim(0, 105)
    coverage_axis.set_xlabel("Teeliste coverage (%)")
    coverage_axis.set_title("Resolved to genus or finer")
    coverage_axis.grid(axis="x", alpha=0.2)

    figure.suptitle(title, fontsize=16)
    figure.tight_layout(rect=[0, 0, 1, 0.96])
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--teeliste", required=True, type=Path)
    parser.add_argument("--db-counts", required=True, type=Path)
    consensus_group = parser.add_mutually_exclusive_group(required=True)
    consensus_group.add_argument("--consensus", type=Path)
    consensus_group.add_argument(
        "--consensus-dir",
        type=Path,
        help="Directory of *.cluster_consensus.tsv files for an all-primer plot",
    )
    parser.add_argument(
        "--db-name",
        help="Database suffix to remove from primer labels (for example refseq_plastid)",
    )
    parser.add_argument("--output", type=Path, default=Path("commnet.html"))
    parser.add_argument("--png", type=Path, help="Optional static PNG output")
    parser.add_argument("--title", default="Teeliste genus-level coverage")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.consensus_dir:
        plot_data = summarize_primers(
            args.teeliste,
            args.db_counts,
            args.consensus_dir,
            args.db_name,
        )
        figure = make_primer_comparison(plot_data, args.title)
    else:
        plot_data = build_plot_data(args.teeliste, args.db_counts, args.consensus)
        figure = make_figure(plot_data, args.title)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(args.output, include_plotlyjs=True)
    data_output = args.output.with_suffix(".tsv")
    plot_data.drop(
        columns=["database_log", "amplified_log"], errors="ignore"
    ).to_csv(data_output, sep="\t", index=False)
    if args.png:
        args.png.parent.mkdir(parents=True, exist_ok=True)
        if args.consensus_dir:
            write_primer_comparison_png(plot_data, args.title, args.png)
        else:
            write_png(plot_data, args.title, args.png)

    print(f"Wrote plot: {args.output}")
    print(f"Wrote plot data: {data_output}")
    if args.png:
        print(f"Wrote image: {args.png}")


if __name__ == "__main__":
    main()
