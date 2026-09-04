# Main Analysis Workflow

The main analysis lives in `workflows/barbeque.nf`. It runs after `DATABASE` has produced the selected database FASTA channel.

## Inputs

Primer input is one of:

- `--primer_set`: named primer sets resolved from the FooDMe2 catalog.
- `--input <samplesheet.tsv>`: a TSV with `primer`, `fwd`, `rev`, `min`, and `max` columns.
- `--input <primers.fasta>` or `--input <directory>`: one `.fa`/`.fasta`/`.fna` primer FASTA, or a directory of them. Requires `--primer_min` and `--primer_max`.

A file input is read as a FASTA or as a samplesheet based on its first non-blank line, not its
extension (`WorkflowPipeline.isFastaInput`). FASTA input - and every `--primer_set` download - is
converted to a samplesheet by `modules/parse_primers` first, so all three routes reach `INPUT_CHECK`
in the same shape and share one validation and staging step.

See [Primer Input](primer_input.md) for the parsing rules and set naming.

Database input is one of:

- `--dbs`: comma-separated installed database ids from `conf/resources.config`.
- `--custom_db`: a user-provided FASTA.

Before primer benchmarking, `DATABASE` can restrict records with `--taxid` and
apply header/length cleaning to every selected database with `--db_filter`.

Taxonomy input is resolved from:

- `--taxdump`, otherwise `<reference_base>/barbeque/<reference_version>/new_taxdump`
- `--accession_taxonomy`, otherwise an accession mapping under `<reference_base>/barbeque/<reference_version>/genbank2taxid`

## Workflow Steps

1. Resolve primers into a common structure: primer id, forward primer, reverse primer, minimum amplicon length, maximum amplicon length.
2. Combine every primer with every selected database.
3. Run in-silico PCR with OBI-PCR.
4. Parse the raw OBI-PCR output.
5. Optionally remove accessions listed by `--accession_blocklist` from both parsed results and amplicon FASTA.
6. Drop primer/database pairs with no retained amplicons from downstream analysis.
7. Optionally apply `--mask` to mimic single-end or paired-end read coverage.
8. Write amplicon length summaries.
9. Build accession-to-taxid tables once per database.
10. Cluster amplicons with `vsearch --cluster_fast` using `--cluster_id`.
11. Parse `.uc` cluster assignments.
12. Join clustered accessions to taxids.
13. Calculate consensus taxonomy per cluster.
14. Summarize database taxonomic distribution.
15. Optionally run target-taxon coverage with `--taxon`.
16. Build one MultiQC report per primer/database combination.

## Excluding Unwanted Accessions

`--accession_blocklist` is applied immediately after OBI-PCR parsing. Matching
accessions are removed from the parsed TSV and amplicon FASTA before any masking,
length profiling, clustering, taxonomy assignment, or reporting. See
[usage.md](usage.md).

## Interactive Results

`--interactive` launches the Streamlit workflow after consensus files are produced. It reads the finished output directory and serves the dashboard on the configured Streamlit port.
