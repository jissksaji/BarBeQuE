#!/usr/bin/env python3
"""
taxid_tree.py

Given a taxid, builds its descendant taxonomic tree using ete3, reports how
many distinct ranks are present, and renders it as an interactive Plotly
sunburst chart. Click a wedge to zoom into that branch, click the center to
zoom back out to the root. No Dash needed - this is native Plotly behavior,
baked into the saved HTML.

Usage:
    python taxid_tree.py --taxid 3398 --output tree.html
    python taxid_tree.py --taxid 4442 --rank_limit species --output camellia.html
"""

import argparse
from ete3 import NCBITaxa
import plotly.graph_objects as go


def parse_args():
    parser = argparse.ArgumentParser(description="Build an interactive taxonomic tree for a given taxid")
    parser.add_argument("--taxid", required=True, type=int, help="Root taxid to build the descendant tree from")
    parser.add_argument("--dbfile", default=None, help="Path to ete3 taxa.sqlite (default: ete3's usual ~/.etetoolkit/taxa.sqlite)")
    parser.add_argument("--rank_limit", default=None, help="Optional rank to stop descending at, e.g. 'species' (skips subspecies/strains)")
    parser.add_argument("--output", default="taxid_tree.html", help="Output HTML file for the interactive chart")
    return parser.parse_args()


def get_descendant_tree(ncbi, taxid, rank_limit):
    """
    Get the full descendant tree under a taxid as an ete3 Tree object.
    collapse_subspecies=False to avoid taxid mismatches (established issue
    with True in this pipeline).
    """
    tree = ncbi.get_descendant_taxa(
        taxid,
        intermediate_nodes=True,
        rank_limit=rank_limit,
        collapse_subspecies=False,
        return_tree=True,
    )
    return tree


def collect_ranks(tree, ncbi):
    """
    Walk every node in the tree and collect the set of distinct ranks present.
    Returns (ranks_present, rank_dict) where rank_dict maps taxid -> rank.
    """
    taxids = [int(node.name) for node in tree.traverse()]
    rank_dict = ncbi.get_rank(taxids)
    ranks_present = set(rank_dict.values())
    return ranks_present, rank_dict


def build_sunburst_data(tree, ncbi, rank_dict):
    """
    Walk the ete3 tree and build the ids, labels, and parents lists needed
    for a Plotly sunburst chart. node.name holds the taxid as a string.
    """
    taxids = [int(node.name) for node in tree.traverse()]
    name_dict = ncbi.get_taxid_translator(taxids)

    ids = []
    labels = []
    parents = []

    for node in tree.traverse():
        taxid = int(node.name)
        rank = rank_dict.get(taxid, "no rank")
        sci_name = name_dict.get(taxid, str(taxid))

        ids.append(node.name)
        labels.append(f"{sci_name} ({rank})")
        parents.append("" if node.is_root() else node.up.name)

    return ids, labels, parents


def make_sunburst(ids, labels, parents, output_path, root_label, n_ranks):
    """
    Build and save an interactive Plotly sunburst chart. Clicking a wedge
    zooms into that branch; clicking the center zooms back out.
    """
    fig = go.Figure(go.Sunburst(
        ids=ids,
        labels=labels,
        parents=parents,
    ))
    fig.update_layout(
        title=f"Taxonomic tree rooted at {root_label} ({n_ranks} distinct ranks)",
        margin=dict(t=40, l=0, r=0, b=0),
    )
    fig.write_html(output_path)
    print(f"Saved interactive tree to {output_path}")


def main():
    args = parse_args()

    ncbi = NCBITaxa(dbfile=args.dbfile) if args.dbfile else NCBITaxa()

    print(f"Building descendant tree for taxid {args.taxid} ...")
    tree = get_descendant_tree(ncbi, args.taxid, args.rank_limit)

    ranks_present, rank_dict = collect_ranks(tree, ncbi)
    print(f"Found {len(ranks_present)} distinct ranks: {sorted(ranks_present)}")
    print(f"Total nodes in tree: {len(list(tree.traverse()))}")

    ids, labels, parents = build_sunburst_data(tree, ncbi, rank_dict)

    root_name = ncbi.get_taxid_translator([args.taxid]).get(args.taxid, str(args.taxid))
    make_sunburst(ids, labels, parents, args.output, root_name, len(ranks_present))


if __name__ == "__main__":
    main()