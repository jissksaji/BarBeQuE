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
    "blocklist_filter"
    "cat_fastq"
    "custom"
    "cutadapt"
    "db_distribution"
    "download"
    "fastp"
    "gunzip"
    "helper"
    "mask"
    "multiqc"
    "obipcr"
    "parse_obipcr"
    "primer_disambiguate"
    "samtools"
    "seqkit"
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

echo "All tests ran successfully!"
