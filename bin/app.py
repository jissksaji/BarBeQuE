import sys
from pathlib import Path
 
import pandas as pd
import plotly.express as px
import streamlit as st
import re
 
st.set_page_config(page_title="Primer Cluster Taxonomy Dashboard", layout="wide")
 
# Usage: python -m streamlit run app.py -- /path/to/results_dir
if len(sys.argv) < 2:
    st.error("Usage: streamlit run app.py -- <path_to_results_dir>")
    st.stop()
 
DATA_DIR = Path(sys.argv[1])
if not DATA_DIR.exists():
    st.error(f"Directory not found: {DATA_DIR}")
    st.stop()

COLUMNS = [
    "cluster_id",
    "accession",
    "accession_taxid",
    "accession_name",
    "assigned_taxid",
    "assigned_rank",
    "assigned_name"
]

DTYPES = {
    "cluster_id": "Int64",
    "accession": "string",
    "accession_taxid": "Int64",
    "accession_name": "string",
    "assigned_taxid": "Int64",
    "assigned_rank": "string",
    "assigned_name": "string"
}


# Helper functions

def load_tsv_files(outdir):
    dfs = []
    consensus_dir = Path(outdir) / "consensus"
    if not consensus_dir.exists():
        st.error(f"Directory not found: {consensus_dir}")
        st.stop()
    for file in sorted(Path(consensus_dir).glob("*.tsv")):
        tmp = pd.read_csv(
            file,
            sep="\t",
            header=None,
            names=COLUMNS,
            dtype=DTYPES
        )
        tmp["primer"] = Path(file.stem).stem

        dfs.append(tmp)

    if not dfs:
        return pd.DataFrame(columns=COLUMNS + ["primer"])

    return pd.concat(dfs, ignore_index=True)


def most_common(series, default="NA"):
    series = series.dropna()
    return series.value_counts().idxmax() if len(series) > 0 else default


def clean_taxon_name(name, rank=None):
    if pd.isna(name):
        return "Unknown"

    name = str(name).split(";")[0].strip()

    if rank == "genus":
        return name.split()[0] if name else "Unknown"

    if rank == "species":
        parts = name.split()
        return " ".join(parts[:2]) if len(parts) >= 2 else name

    return name


def shorten(text, max_len=35):
    text = str(text)
    return text if len(text) <= max_len else text[:max_len] + "..."


def compact_plot(fig, height=360):
    fig.update_layout(
        height=height,
        margin=dict(l=20, r=20, t=60, b=90),
        xaxis_tickangle=-45
    )
    return fig


# Load data

st.title("Primer Cluster Taxonomy Dashboard")

df = load_tsv_files(DATA_DIR)

if df.empty:
    st.warning("No TSV files found in the selected folder.")
    st.stop()

# Primer filter
all_primers = sorted(df["primer"].unique())
primer_choice = st.sidebar.selectbox("Filter by primer", ["All"] + all_primers)

if primer_choice != "All":
    df = df[df["primer"] == primer_choice]


# Basic statistics
st.subheader("Basic statistics")

total_rows = len(df)
total_clusters = df["cluster_id"].nunique()
total_accessions = df["accession"].nunique()
total_accession_taxids = df["accession_taxid"].nunique()
total_assigned_taxids = df["assigned_taxid"].nunique()

most_abundant_taxid = most_common(df["assigned_taxid"])
most_common_rank = most_common(df["assigned_rank"])

matching_rows = df.loc[
    df["assigned_taxid"] == most_abundant_taxid,
    ["assigned_name", "assigned_rank"]
].dropna()

if len(matching_rows) > 0:
    row = matching_rows.iloc[0]
    most_abundant_name = clean_taxon_name(
        row["assigned_name"],
        row["assigned_rank"]
    )
else:
    most_abundant_name = "Unknown"

col1, col2, col3, col4 = st.columns(4)

col1.metric("Total clusters", total_clusters)

col2.metric("Total accessions", total_accessions)
col2.metric("Unique accession taxids", total_accession_taxids)

col3.metric("Unique assigned taxids", total_assigned_taxids)
col3.metric("Most common assigned rank", most_common_rank)

col4.metric(
    "Most abundant assigned taxon",
    most_abundant_name,
    delta=f"taxid {most_abundant_taxid} · {most_common(df.loc[df['assigned_taxid'] == most_abundant_taxid, 'assigned_rank'])}",
    delta_color="off"
)


# Cluster summary
st.subheader("Sequences per cluster")

seqs_per_cluster = (
    df.groupby(["primer", "cluster_id"], dropna=False)
    .agg(
        sequence_count=("accession", "count"),
        unique_accession_taxids=("accession_taxid", "nunique"),
        assigned_taxid=("assigned_taxid", "first"),
        assigned_rank=("assigned_rank", "first"),
        assigned_name=("assigned_name", "first")
    )
    .reset_index()
    .sort_values("sequence_count", ascending=False)
)

st.dataframe(seqs_per_cluster, width="stretch", height=300)


# Rank and primer summaries
left, right = st.columns(2)

with left:
    st.subheader("Assigned rank distribution")

    rank_counts = (
        seqs_per_cluster["assigned_rank"]
        .fillna("unresolved")   # ← cluster level, consistent
        .value_counts()
        .reset_index()
    )
    rank_counts.columns = ["assigned_rank", "count"]

    st.dataframe(rank_counts, width="stretch", height=220)

with right:
    st.subheader("Primer-level summary")

    primer_summary = (
        df.groupby("primer")
        .agg(
            total_rows=("accession", "count"),
            total_clusters=("cluster_id", "nunique"),
            total_accessions=("accession", "nunique"),
            unique_accession_taxids=("accession_taxid", "nunique"),
            unique_assigned_taxids=("assigned_taxid", "nunique"),
            most_common_rank=(
                "assigned_rank",
                lambda x: x.dropna().mode().iloc[0]
                if not x.dropna().mode().empty
                else "NA"
            )
        )
        .reset_index()
    )

    st.dataframe(primer_summary, width="stretch", height=220)


# Most abundant assigned taxa
st.subheader("Most abundant assigned taxa")

top_n = st.slider("Number of taxa to show", 5, 50, 20)

taxon_counts = (
    seqs_per_cluster
    .groupby(["assigned_taxid", "assigned_rank"], dropna=False)
    .agg(
        count=("cluster_id", "count"),
        assigned_name=("assigned_name", "first")
    )
    .reset_index()
    .sort_values("count", ascending=False)
)

taxon_counts["clean_name"] = taxon_counts.apply(
    lambda row: clean_taxon_name(row["assigned_name"], row["assigned_rank"]),
    axis=1
)

taxon_counts["plot_name"] = taxon_counts.apply(
    lambda row: f'{row["assigned_taxid"]} | {shorten(row["clean_name"], 30)}',
    axis=1
)

st.dataframe(taxon_counts, width="stretch", height=300)

top_taxa = taxon_counts.head(top_n)

fig_taxa = px.bar(
    top_taxa,
    x="plot_name",
    y="count",
    color="assigned_rank",
    title=f"Top {top_n} most abundant assigned taxa by cluster count",
    hover_data=[
        "assigned_taxid",
        "assigned_rank",
        "clean_name",
        "count"
    ],
    category_orders={"plot_name": top_taxa["plot_name"].tolist()}
)

st.plotly_chart(
    compact_plot(fig_taxa, height=450),
    width="stretch"
)


# Cluster resolution viewer
st.subheader("Cluster resolution viewer")

viewer_primers = sorted(df["primer"].unique())
viewer_primer = st.selectbox("Select a primer", viewer_primers)
viewer_df = df[df["primer"] == viewer_primer]

selected_cluster = st.selectbox(
    "Select a cluster",
    sorted(viewer_df["cluster_id"].dropna().unique())
)

cluster_df = viewer_df[viewer_df["cluster_id"] == selected_cluster]

cluster_rank = most_common(cluster_df["assigned_rank"])
cluster_taxid = most_common(cluster_df["assigned_taxid"])

cluster_name_rows = cluster_df.loc[
    cluster_df["assigned_taxid"] == cluster_taxid,
    ["assigned_name", "assigned_rank"]
].dropna()

if len(cluster_name_rows) > 0:
    cluster_name = clean_taxon_name(
        cluster_name_rows.iloc[0]["assigned_name"],
        cluster_name_rows.iloc[0]["assigned_rank"]
    )
else:
    cluster_name = "Unknown"

c1, c2, c3 = st.columns(3)
c1.metric("Selected cluster", selected_cluster)
c2.metric("Resolved taxid", cluster_taxid)
c3.metric("Resolved rank", cluster_rank)

st.success(f"Resolved taxon: {cluster_name}")

taxids_in_cluster = (
    cluster_df
    .groupby(["accession_taxid", "accession_name"], dropna=False)
    .size()
    .reset_index(name="sequence_count")
    .sort_values("sequence_count", ascending=False)
)

taxids_in_cluster["plot_name"] = taxids_in_cluster.apply(
    lambda row: f'{row["accession_taxid"]} | {shorten(row["accession_name"], 30)}',
    axis=1
)

st.dataframe(cluster_df, width="stretch", height=250)

fig_cluster_taxids = px.bar(
    taxids_in_cluster,
    x="plot_name",
    y="sequence_count",
    title=f"Taxids inside cluster {selected_cluster}",
    hover_data=[
        "accession_taxid",
        "accession_name",
        "sequence_count"
    ],
    category_orders={"plot_name": taxids_in_cluster["plot_name"].tolist()}
)

st.plotly_chart(
    compact_plot(fig_cluster_taxids, height=400),
    width="stretch"
)


# Primer resolution efficiency

st.subheader("Primer resolution efficiency (LCA per cluster)")

resolution = (
    seqs_per_cluster
    .groupby(["primer", "assigned_rank"], dropna=False)
    .agg(cluster_count=("cluster_id", "count"))
    .reset_index()
)

# total clusters per primer
primer_totals = (
    seqs_per_cluster
    .groupby("primer")
    .agg(total_clusters=("cluster_id", "count"))
    .reset_index()
)

resolution = resolution.merge(primer_totals, on="primer")
resolution["pct"] = (
    resolution["cluster_count"] / resolution["total_clusters"] * 100
).round(1)

resolution["assigned_rank"] = resolution["assigned_rank"].fillna("unresolved")


st.markdown("**Total clusters per primer**")
cols = st.columns(len(primer_totals))
for col, (_, row) in zip(cols, primer_totals.iterrows()):
    col.metric(row["primer"], f'{int(row["total_clusters"]):,} clusters')
fig_res = px.bar(
    resolution,
    x="primer",
    y="pct",
    color="assigned_rank",
    text="cluster_count",
    title="% of clusters resolved per LCA rank per primer",
    labels={"pct": "% of clusters", "cluster_count": "# clusters"},
    barmode="stack"
)

fig_res.update_traces(textposition="inside")

st.plotly_chart(compact_plot(fig_res, height=450), width="stretch")



#taxonomic coverage from module read directly
def load_tax_coverage(data_dir):
    dfs = []
    tax_dir = Path(data_dir) / "tax_coverage"

    for file in sorted(tax_dir.glob("*.tax_coverage.tsv")):
        tmp = pd.read_csv(
            file,
            sep="\t",
            header=None,
            names=["Taxon", "Status", "Taxid", "Color"],
            skiprows=1
        )
        stem  = file.stem.replace(".tax_coverage", "")
        parts = stem.split("--")
        tmp["primer"] = parts[0].rstrip("_")
        tmp["taxon"]  = parts[1] if len(parts) > 1 else "Unknown"
        dfs.append(tmp)

    if not dfs:
        return pd.DataFrame(columns=["Taxon", "Status", "Taxid", "Color", "primer"])

    return pd.concat(dfs, ignore_index=True)



st.divider()

tax_df = load_tax_coverage(DATA_DIR)

if tax_df.empty:
    st.info("No taxonomic coverage data found. Run pipeline with --taxon to enable.")

else:
    db_species = tax_df[tax_df["Status"] != "NO_DATA"]
    taxon_name = tax_df["taxon"].iloc[0]
    st.header(f"Taxonomic Coverage — {taxon_name}")


    # coverage per primer
    coverage_rows = []
    for primer, group in db_species.groupby("primer"):
        total     = len(group)
        recovered = len(group[group["Status"] == "OK"])
        coverage  = round(recovered / total * 100, 1) if total > 0 else 0.0
        coverage_rows.append({
            "primer":    primer,
                "total":     total,
                "recovered": recovered,
                "coverage":  coverage
            })
        coverage_df = pd.DataFrame(coverage_rows)

    # metrics row
    cols = st.columns(len(coverage_df))
    for col, (_, row) in zip(cols, coverage_df.iterrows()):
        col.metric(
            row["primer"],
            f'{row["coverage"]}%',
            delta=f'{int(row["recovered"])}/{int(row["total"])} species',
            delta_color="off"
        )

    # coverage bar chart
    fig_cov = px.bar(
        coverage_df,
        x="primer",
        y="coverage",
        color="primer",
        title="Taxonomic Coverage per Primer (%)",
        labels={"coverage": "Coverage (%)"}
    )
    st.plotly_chart(compact_plot(fig_cov), use_container_width=True)

    # status breakdown
    status_counts = tax_df.groupby(["primer", "Status"]).size().reset_index(name="count")
    fig_status = px.bar(
        status_counts,
        x="primer",
        y="count",
        color="Status",
        barmode="stack",
        color_discrete_map={
            "OK": "#7ee076",
            "FAIL": "#ff4500",
            "NO_DATA": "#eeeeee"
        },
        title="Species Status Breakdown per Primer"
    )
    st.plotly_chart(compact_plot(fig_status), use_container_width=True)