#!/usr/bin/env bash
set -euo pipefail

cluster_accessions="$1"
accession_taxid="$2"
out_taxids="$3"
out_missing="$4"

awk -F'\t' '
    BEGIN {
        OFS = "\t"
    }

    # First file: cluster_accessions.tsv
    # cluster_id accession
    NR == FNR {
        cluster[$2] = $1
        wanted[$2] = 1
        next
    }

    # Second file: accession_taxid.tsv
    # accession taxid

    # Skip header if present
    FNR == 1 && tolower($1) == "accession" {
        next
    }

    $1 in cluster {
        acc = $1
        taxid = $2
        cid = cluster[acc]

        found[acc] = 1

        key = cid SUBSEP taxid

        # avoid duplicate taxid inside same cluster
        if (!(key in seen)) {
            seen[key] = 1

            if (cid in taxids) {
                taxids[cid] = taxids[cid] " " taxid
            } else {
                taxids[cid] = taxid
            }
        }
    }

    END {
        for (cid in taxids) {
            print cid, taxids[cid] > "'"$out_taxids"'"
        }

        for (acc in wanted) {
            if (!(acc in found)) {
                print cluster[acc], acc > "'"$out_missing"'"
            }
        }
    }
' "$cluster_accessions" "$accession_taxid"

# Ensure files exist even if empty
touch "$out_taxids"
touch "$out_missing"

sort -k1,1n "$out_taxids" -o "$out_taxids"
sort -k1,1n -k2,2 "$out_missing" -o "$out_missing"