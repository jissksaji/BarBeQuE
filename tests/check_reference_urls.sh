#!/usr/bin/env bash
#
# Fast smoke test: verifies every URL declared in conf/reference_sources.config
# is reachable (HTTP HEAD only - no actual download). Safe to run in CI to catch
# link rot (moved/renamed database releases) without downloading anything.
#
# Usage: tests/check_reference_urls.sh
set -euo pipefail

cd "$(dirname "$0")/.."

fail=0

while IFS= read -r line; do
    id=$(echo "$line" | sed -E "s/^params\.reference_sources\.([^.]+)\.url.*/\1/")
    url=$(echo "$line" | sed -E "s/^params\.reference_sources\.[^.]+\.url = '(.*)'$/\1/")
    status=$(curl -o /dev/null -s -L --max-time 20 -w '%{http_code}' -I "$url" || echo "000")
    if [[ "$status" =~ ^(200|301|302)$ ]]; then
        printf 'OK    %-32s %s (%s)\n' "$id" "$url" "$status"
    else
        printf 'FAIL  %-32s %s (%s)\n' "$id" "$url" "$status"
        fail=1
    fi
done < <(nextflow config -flat . | grep -E '^params\.reference_sources\.[^.]+\.url =')

echo
if [[ $fail -eq 0 ]]; then
    echo "All reference_sources URLs are reachable."
else
    echo "One or more reference_sources URLs FAILED - see above."
fi
exit $fail
