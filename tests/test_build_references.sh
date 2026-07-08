#!/usr/bin/env bash
#
# End-to-end test: actually runs --build_references into a scratch directory and
# verifies that every database wired up for automated download (i.e. every entry
# in conf/resources.config that declares a `source`) produces a non-empty file at
# the exact path declared as its `db` value. Pass --install-taxdump / --install-genbank
# to also exercise those two (large, normally opt-in) downloads.
#
# This does real downloads and can take a long time depending on bandwidth.
#
# Usage: tests/test_build_references.sh [--install-taxdump] [--install-genbank]
set -euo pipefail

cd "$(dirname "$0")/.."

ref_base=$(mktemp -d)
trap 'rm -rf "$ref_base"' EXIT

install_taxdump=false
install_genbank=false
extra_args=()
for arg in "$@"; do
    case "$arg" in
        --install-taxdump) extra_args+=(--install_taxdump); install_taxdump=true ;;
        --install-genbank) extra_args+=(--install_genbank); install_genbank=true ;;
        *) echo "Unknown argument: $arg" >&2; exit 2 ;;
    esac
done

echo "Building references into: $ref_base"
nextflow run main.nf --build_references --reference_base "$ref_base" "${extra_args[@]}"

# `nextflow config` does not accept ad hoc --param overrides (only `nextflow run` does),
# so we resolve paths with reference_base unset (-> literal "null/..." prefix) and swap
# in our scratch dir ourselves, rather than re-invoking config with --reference_base.
cfg=$(nextflow config -flat .)

resolve_path() {
    echo "$1" | sed -E "s|^null/|${ref_base}/|"
}

fail=0

while IFS= read -r db_id; do
    path=$(resolve_path "$(echo "$cfg" | grep "^params.references.databases.${db_id}.db " | sed -E "s/.*= '(.*)'\$/\1/")")
    if [[ -s "$path" ]]; then
        printf 'OK    %-20s %s\n' "$db_id" "$path"
    else
        printf 'FAIL  %-20s %s (missing or empty)\n' "$db_id" "$path"
        fail=1
    fi
done < <(echo "$cfg" | grep -E '^params\.references\.databases\.[^.]+\.source =' | sed -E 's/^params\.references\.databases\.([^.]+)\.source.*/\1/')

taxdump_path=$(resolve_path "$(echo "$cfg" | grep '^params.references.taxdump ' | sed -E "s/.*= '(.*)'\$/\1/")")

if $install_taxdump; then
    if [[ -f "${taxdump_path}/nodes.dmp" ]]; then
        printf 'OK    %-20s %s\n' "taxdump" "$taxdump_path"
    else
        printf 'FAIL  %-20s %s (nodes.dmp missing)\n' "taxdump" "$taxdump_path"
        fail=1
    fi
fi

if $install_genbank; then
    genbank_path="${taxdump_path}/nucl_gb.accession2taxid"
    if [[ -s "$genbank_path" ]]; then
        printf 'OK    %-20s %s\n' "genbank" "$genbank_path"
    else
        printf 'FAIL  %-20s %s (missing or empty)\n' "genbank" "$genbank_path"
        fail=1
    fi
fi

echo
if [[ $fail -eq 0 ]]; then
    echo "All wired-up databases downloaded successfully."
else
    echo "One or more databases FAILED to download - see above."
fi
exit $fail
