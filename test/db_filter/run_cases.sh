#!/bin/bash
# Exercise every combination of the DB_FILTER parameters against test/db_filter/reference.fasta.
#
# Each filter is independently optional, so each must be checked on its own as well as
# combined - a shared gate silently disabled two of them once already.
#
# Usage: test/db_filter/run_cases.sh [-profile conda]

set -u

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
profile="${*:--profile conda}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

failures=0

# name | min_length | max_length | max_n | expected surviving record IDs
cases=(
  "no optional filters|null|null|null|clean,short20,long400,amb10,amb11,lower_amb11,homopolymer"
  "min length only|50|null|null|clean,long400,amb10,amb11,lower_amb11,homopolymer"
  "max length only|null|100|null|clean,short20,amb10,amb11,lower_amb11,homopolymer"
  "both length bounds|50|100|null|clean,amb10,amb11,lower_amb11,homopolymer"
  "ambiguity only|null|null|10|clean,short20,long400,amb10,homopolymer"
  "pipeline defaults|50|null|10|clean,long400,amb10,homopolymer"
  "strictest ambiguity|null|null|0|clean,short20,long400,homopolymer"
)

for entry in "${cases[@]}"; do
  IFS='|' read -r name min max maxn expect <<< "$entry"
  cat > "$work/params.json" <<JSON
{
  "db_filter_min_length": $min,
  "db_filter_max_length": $max,
  "db_filter_max_n": $maxn,
  "expect": "$expect"
}
JSON
  if nextflow run "$here/main.nf" $profile -params-file "$work/params.json" \
        -w "$work/w" > "$work/log.txt" 2>&1; then
    printf 'PASS  %s\n' "$name"
  else
    printf 'FAIL  %s\n' "$name"
    sed -n '/DB_FILTER kept the wrong records/,/wrongly kept/p' "$work/log.txt" | sed 's/^/        /'
    failures=$((failures + 1))
  fi
  rm -rf "$work/w"
done

if [ "$failures" -ne 0 ]; then
  echo "$failures of ${#cases[@]} DB_FILTER cases failed."
  exit 1
fi

echo "All ${#cases[@]} DB_FILTER cases passed."
