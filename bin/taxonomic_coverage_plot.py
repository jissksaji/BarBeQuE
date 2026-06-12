#!/usr/bin/env python3
"""
taxonomic_coverage_plot.py — BarBeQuE Taxonomic Coverage Plot

Usage:
    taxonomic_coverage_plot.py --input <tax_coverage.tsv> [<...>] --output <out.html>
"""

import argparse
import warnings
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from pathlib import Path


STATUS_COLORS = {"OK": "#4daf4a", "FAIL": "#e41a1c", "NO_DATA": "#bdbdbd"}
STATUS_ORDER  = ["OK", "FAIL", "NO_DATA"]
RANK_ORDER    = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "no rank"]

COMMON_AXIS = dict(
    showline=True, linewidth=2, linecolor="black",
    ticks="outside", tickwidth=2, tickcolor="black", ticklen=6,
    tickfont=dict(size=13, family="Arial", color="black"),
    title_font=dict(size=15, family="Arial", color="black"),
)
Y_GRID = dict(gridcolor="rgba(200,200,200,0.4)", gridwidth=1)

LEGEND_STYLE = dict(
    title="<b>Status</b>",
    orientation="h",
    yanchor="bottom", y=1.04, xanchor="right", x=1,
    font=dict(size=13, family="Arial", color="black"),
    bgcolor="rgba(255,255,255,0.9)",
    bordercolor="black", borderwidth=1.5,
)


# ── I/O ──────────────────────────────────────────────────────────────────────

def read_tax_coverage(paths):
    dfs = []
    for path in paths:
        p     = Path(path)
        stem  = p.name.replace(".tsv", "").replace(".tax_coverage", "")
        parts = stem.split("--")

        if len(parts) < 2:
            warnings.warn(
                f"Cannot parse primer/taxon from '{p.name}'; "
                f"expected '<primer>--<taxon>.tax_coverage.tsv'. "
                f"Using full stem as primer, 'Unknown' as taxon."
            )

        primer     = parts[0].rstrip("_") if parts else stem
        taxon_name = parts[1] if len(parts) > 1 else "Unknown"

        # Auto-detect header: skip row 0 if first cell looks like a column label
        raw        = pd.read_csv(p, sep="\t", header=None, nrows=1)
        first_cell = str(raw.iloc[0, 0]).strip().lower()
        skip       = 1 if first_cell in ("taxon", "name", "species", "#taxon") else 0

        df = pd.read_csv(
            p, sep="\t", header=None, skiprows=skip,
            names=["Taxon", "Status", "Taxid", "Color", "Rank"]
        )
        df["Rank"]       = df["Rank"].fillna("species") if "Rank" in df.columns else "species"
        df["primer"]     = primer
        df["taxon_name"] = taxon_name
        dfs.append(df)

    if not dfs:
        return pd.DataFrame(
            columns=["Taxon", "Status", "Taxid", "Color", "Rank", "primer", "taxon_name"]
        )
    return pd.concat(dfs, ignore_index=True)


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_sorted_primers(status_counts):
    """Primers sorted by OK% descending so the best primer is always first."""
    totals = status_counts.groupby("primer")["count"].sum()
    ok_cnt = status_counts[status_counts["Status"] == "OK"].set_index("primer")["count"]
    ok_pct = (ok_cnt.reindex(totals.index, fill_value=0) / totals * 100).sort_values(ascending=False)
    return ok_pct.index.tolist()


def _trunc(label, n=30):
    return label if len(label) <= n else label[:13] + "…" + label[-(n - 14):]


# ── Figure 1: absolute + normalized side-by-side ──────────────────────────────

def build_coverage_figure(tax_df):
    status_counts = (
        tax_df.groupby(["primer", "Status"])
              .size()
              .reset_index(name="count")
    )
    totals = status_counts.groupby("primer")["count"].transform("sum")
    status_counts["pct"] = (status_counts["count"] / totals * 100).round(1)
    status_counts["n"]   = totals

    sorted_primers = get_sorted_primers(status_counts)
    taxon_str      = ", ".join(tax_df["taxon_name"].unique()) or "Unknown"

    # Lookup tables indexed by primer for fast access
    n_per_primer = (
        status_counts.groupby("primer")["n"].first()
    )
    max_n = int(n_per_primer.max())

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["<b>Absolute counts</b>", "<b>Relative (%)</b>"],
        horizontal_spacing=0.12,
    )

    for status in STATUS_ORDER:
        sub    = status_counts[status_counts["Status"] == status].set_index("primer")
        ys_abs = [int(sub.loc[p, "count"]) if p in sub.index else 0   for p in sorted_primers]
        ys_pct = [float(sub.loc[p, "pct"]) if p in sub.index else 0.0 for p in sorted_primers]
        ns     = [int(n_per_primer[p]) for p in sorted_primers]

        # Only render in-bar text when the segment is wide enough to be readable
        text_abs = [
            f"<b>{c}</b><br><span style='font-size:11px'>({pct:.0f}%)</span>"
            if n > 0 and c / n > 0.05 else ""
            for c, n, pct in zip(ys_abs, ns, ys_pct)
        ]
        text_pct = [
            f"<b>{pct:.0f}%</b>" if pct > 5 else ""
            for pct in ys_pct
        ]

        shared_bar = dict(
            marker_color=STATUS_COLORS[status],
            marker_line_color="black",
            marker_line_width=1.2,
            textposition="inside",
            insidetextanchor="middle",
            legendgroup=status,
        )

        # Absolute
        fig.add_trace(
            go.Bar(
                name=status,
                x=sorted_primers, y=ys_abs,
                text=text_abs,
                customdata=list(zip(ys_abs, ns, ys_pct)),
                hovertemplate=(
                    "<b>%{x}</b><br>"
                    f"Status: {status}<br>"
                    "Count: %{customdata[0]} / %{customdata[1]}<br>"
                    "Percentage: %{customdata[2]:.1f}%<extra></extra>"
                ),
                showlegend=True,
                **shared_bar,
            ),
            row=1, col=1,
        )

        # Normalized
        fig.add_trace(
            go.Bar(
                name=status,
                x=sorted_primers, y=ys_pct,
                text=text_pct,
                customdata=list(zip(ys_abs, ns, ys_pct)),
                hovertemplate=(
                    "<b>%{x}</b><br>"
                    f"Status: {status}<br>"
                    "%{customdata[2]:.1f}% (%{customdata[0]} / %{customdata[1]})<extra></extra>"
                ),
                showlegend=False,
                **shared_bar,
            ),
            row=1, col=2,
        )

    # n= annotation above each absolute bar
    for primer in sorted_primers:
        n = int(n_per_primer[primer])
        fig.add_annotation(
            x=primer, y=n + max_n * 0.03,
            text=f"<i>n={n}</i>",
            showarrow=False,
            font=dict(size=11, family="Arial", color="#555"),
            row=1, col=1,
        )

    tick_texts = [_trunc(p) for p in sorted_primers]
    tick_angle = 0 if len(sorted_primers) == 1 else -45

    fig.update_layout(
        barmode="stack",
        title=dict(
            text=(
                "<b>Species Amplification Status per Primer</b>"
                f"<br><span style='font-size:14px; color:#555'>{taxon_str} · sorted by OK%</span>"
            ),
            y=0.97, x=0.5, xanchor="center", yanchor="top",
            font=dict(size=20, family="Arial", color="black"),
        ),
        legend=LEGEND_STYLE,
        font=dict(family="Arial", size=13, color="black"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=700,
        margin=dict(t=140, b=170, l=80, r=40),
    )

    for col in (1, 2):
        fig.update_xaxes(
            **COMMON_AXIS,
            tickmode="array",
            tickvals=sorted_primers,
            ticktext=tick_texts,
            tickangle=tick_angle,
            row=1, col=col,
        )
    fig.update_yaxes(**COMMON_AXIS, **Y_GRID, title_text="<b>Number of Species</b>", row=1, col=1)
    fig.update_yaxes(**COMMON_AXIS, **Y_GRID, title_text="<b>Species (%)</b>",        row=1, col=2)

    return fig, sorted_primers


# ── Figure 2: status by resolved rank ─────────────────────────────────────────

def build_rank_figure(tax_df, sorted_primers):
    """One subplot per primer showing Status × Rank.

    Distinguishes 'primer failed to amplify' from 'amplified but classified
    only to genus/family' — directly relevant to rbcL/ITS2 resolution issues.
    Skipped if all entries have the same rank (no information gain).
    """
    if set(tax_df["Rank"].dropna().unique()) <= {"species"}:
        return None

    rank_counts = (
        tax_df.groupby(["primer", "Rank", "Status"])
              .size()
              .reset_index(name="count")
    )

    n_primers = len(sorted_primers)
    fig = make_subplots(
        rows=1, cols=n_primers,
        subplot_titles=[f"<b>{_trunc(p, 25)}</b>" for p in sorted_primers],
        shared_yaxes=True,
        horizontal_spacing=max(0.03, 0.15 / n_primers),
    )

    for col_i, primer in enumerate(sorted_primers, start=1):
        sub           = rank_counts[rank_counts["primer"] == primer]
        present_ranks = [r for r in RANK_ORDER if r in sub["Rank"].values]

        for status in STATUS_ORDER:
            ys = []
            for rank in present_ranks:
                rows = sub[(sub["Rank"] == rank) & (sub["Status"] == status)]
                ys.append(int(rows["count"].iloc[0]) if len(rows) else 0)

            fig.add_trace(
                go.Bar(
                    name=status,
                    x=present_ranks, y=ys,
                    marker_color=STATUS_COLORS[status],
                    marker_line_color="black",
                    marker_line_width=1.0,
                    showlegend=(col_i == 1),
                    legendgroup=status,
                    hovertemplate=(
                        f"<b>{_trunc(primer, 40)}</b><br>"
                        "Rank: %{x}<br>"
                        f"Status: {status}<br>"
                        "Count: %{y}<extra></extra>"
                    ),
                ),
                row=1, col=col_i,
            )

    fig.update_layout(
        barmode="stack",
        title=dict(
            text=(
                "<b>Status by Resolved Rank per Primer</b>"
                "<br><span style='font-size:13px; color:#555'>"
                "Distinguishes amplification failure (FAIL at species) from resolution loss (FAIL at genus/family)"
                "</span>"
            ),
            y=0.97, x=0.5, xanchor="center", yanchor="top",
            font=dict(size=18, family="Arial", color="black"),
        ),
        legend=LEGEND_STYLE,
        font=dict(family="Arial", size=12, color="black"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=500,
        margin=dict(t=130, b=110, l=80, r=40),
    )
    fig.update_xaxes(**COMMON_AXIS, tickangle=-30)
    fig.update_yaxes(**COMMON_AXIS, **Y_GRID, title_text="<b>Count</b>", row=1, col=1)

    return fig


# ── HTML output ───────────────────────────────────────────────────────────────

def write_html(figs, path):
    html_parts = [
        pio.to_html(fig, full_html=False, include_plotlyjs=(i == 0))
        for i, fig in enumerate(figs)
    ]

    full = (
        "<!DOCTYPE html>\n<html>\n<head>\n"
        '  <meta charset="utf-8">\n'
        "  <title>BarBeQuE — Taxonomic Coverage</title>\n"
        "  <style>\n"
        "    body { font-family: Arial, sans-serif; background: white; margin: 0; padding: 20px; }\n"
        "    .section { margin-bottom: 40px; border-top: 2px solid #e8e8e8; padding-top: 20px; }\n"
        "    .section:first-child { border-top: none; }\n"
        "  </style>\n"
        "</head>\n<body>\n"
        + "".join(f'<div class="section">{h}</div>\n' for h in html_parts)
        + "</body>\n</html>"
    )
    Path(path).write_text(full, encoding="utf-8")
    print(f"[taxonomic_coverage_plot] HTML -> {path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="BarBeQuE taxonomic coverage plot")
    parser.add_argument("--input",  nargs="+", required=True, help=".tax_coverage.tsv files")
    parser.add_argument("--output", required=True,            help="Output HTML path")
    args = parser.parse_args()

    df = read_tax_coverage(args.input)
    if df.empty:
        print("[taxonomic_coverage_plot] No data found.")
        return

    coverage_fig, sorted_primers = build_coverage_figure(df)
    figs = [coverage_fig]

    rank_fig = build_rank_figure(df, sorted_primers)
    if rank_fig is not None:
        figs.append(rank_fig)

    write_html(figs, args.output)


if __name__ == "__main__":
    main()