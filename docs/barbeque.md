# Main Analysis Workflow

The main analysis lives in `workflows/barbeque.nf`. It runs after `DATABASE` has produced the selected database FASTA channel.

## Inputs

Primer input is one of:

- `--primer_set`: named primer sets resolved from the FooDMe2 catalog.
- `--input <samplesheet.tsv>`: a TSV with `primer`, `fwd`, `rev`, `min`, and `max` columns.
- `--input <directory>`: `.fa`, `.fasta`, or `.fna` primer FASTAs. Directory mode requires `--primer_min` and `--primer_max`.

Database input is one of:

- `--dbs`: comma-separated installed database ids from `conf/resources.config`.
- `--custom_db`: a user-provided FASTA.

Before primer benchmarking, `DATABASE` can restrict records with `--taxid`,
remove records listed by FooDMe2 with `--blocklist`, and apply header/length
cleaning with `--custom_db_filter`.

Taxonomy input is resolved from:

- `--taxdump`, otherwise `<reference_base>/barbeque/<reference_version>/new_taxdump`
- `--accession_taxonomy`, otherwise an accession mapping under `<reference_base>/barbeque/<reference_version>/genbank2taxid`

## Workflow Steps

1. Resolve primers into a common structure: primer id, forward primer, reverse primer, minimum amplicon length, maximum amplicon length.
2. Combine every primer with every selected database.
3. Run in-silico PCR with `--insilico_tool obipcr` or `--insilico_tool cutadapt`.
4. Parse raw OBI output when `obipcr` is used.
5. Drop primer/database pairs with no amplicons from downstream analysis.
6. Optionally apply `--mask` to mimic single-end or paired-end read coverage.
7. Write amplicon length summaries.
8. Build accession-to-taxid tables once per database.
9. Optionally run advisory species-divergence screening with `--screen_species_divergence`.
10. Optionally remove curated bad accessions with `--exclude_accessions`.
11. Cluster amplicons with `vsearch --cluster_fast` using `--cluster_id`.
12. Parse `.uc` cluster assignments.
13. Join clustered accessions to taxids.
14. Calculate consensus taxonomy per cluster.
15. Summarize database taxonomic distribution.
16. Optionally write completeness tables with `--completeness_table`.
17. Optionally run target-taxon coverage with `--taxon`.
18. Build one MultiQC report per primer/database combination.

## Advisory Screening vs Filtering

`--screen_species_divergence` reports candidate misannotations. It does not change clustering input.

`--exclude_accessions <file>` changes clustering input. Use this after reviewing divergence reports. The file is plain text, one accession per line, with optional `#` comments.

## Interactive Results

`--interactive` launches the Streamlit workflow after consensus files are produced. It reads the finished output directory and serves the dashboard on the configured Streamlit port.
