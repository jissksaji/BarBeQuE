#!/bin/bash
set -e
export PATH="$PWD/bin:$PATH"

echo "Running tests for all modules..."

echo "==================================="
echo "Running Python utility unit tests"
echo "==================================="
python3 -m unittest discover -s tests -p 'test*.py' -v

modules=(
    "blast"
    "cat_fastq"
    "custom"
    "db_distribution"
    "download"
    "fastp"
    "gunzip"
    "helper"
    "mask"
    "multiqc"
    "obipcr"
    "filter_obipcr"
    "parse_obipcr"
    "parse_primers"
    "primer_disambiguate"
    "samtools"
    "streamlit"
    "taxonkit"
    "untar"
    "vsearch"
)

for mod in "${modules[@]}"; do
    echo "==================================="
    echo "Running test for module: $mod"
    echo "==================================="
    nextflow run "test/${mod}/main.nf" -profile conda
done

echo "==================================="
echo "Running DB_FILTER parameter combinations"
echo "==================================="
# Each db_filter param is independently optional, so every combination is checked -
# a single shared gate silently disabled two of the three filters once already.
test/db_filter/run_cases.sh -profile conda

echo "All tests ran successfully!"
