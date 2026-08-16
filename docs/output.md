# Outputs

Outputs are written under `--outdir`, default `results`.

## Core Directories

| Directory | Contents |
| --- | --- |
| `pipeline_info/` | Staged samplesheet and software version metadata. |
| `reports/` | MultiQC reports, one per primer/database group. |
| `raw/obipcr/` | Raw OBI in-silico PCR output when `--insilico_tool obipcr`. |
| `raw/cutadapt/` | Raw cutadapt in-silico PCR output when `--insilico_tool cutadapt`. |
| `parsed_obipcr/` | Parsed OBI results with primer mismatch and amplicon metrics. |
| `collapsed_primers/` | Consensus primer FASTAs created from primer FASTA directory input. |
| `amplicon_lengths/` | Amplicon length summaries. |
| `build_db_taxids/` | Accession-to-taxid and taxid-count tables built from each database. |
| `db_distribution/` | Taxonomic composition summaries for each database. |
| `cluster_fast/` | VSEARCH clustered FASTA and `.uc` files. |
| `join_accession_taxonomy/` | Cluster accession assignments joined to taxids. |
| `consensus/` | Consensus taxonomy per cluster. |

## Optional Directories

| Directory | Created When | Contents |
| --- | --- | --- |
| `taxid_filtered/` | `--taxid` | Database FASTAs restricted to the requested taxon. |
| `blocklist_filtered/` | `--blocklist` | FASTAs with FooDMe2-blocklisted taxids removed, plus TSV removal summaries. |
| `filtered_amplicons/` | `--exclude_accessions` | Amplicon FASTAs after curated accession removal and exclusion audit tables. |
| `species_split/` | `--screen_species_divergence` | Amplicons split into species-level FASTAs. |
| `species_divergence/` | `--screen_species_divergence` | Divergence tables, divergent FASTAs, and linkage matrices. |
| `tax_coverage/` | `--taxon` | Taxon-focused coverage and species representation tables. |
| `completeness/` | `--completeness_table` | Database completeness summaries for `--taxon`. |
| `hierarchical_clustering/` | `--hierarchical_clustering` | Standalone clustering outputs for a custom FASTA. |

## Consensus Table

Files in `consensus/` are the main result tables. They report cluster-level consensus taxonomy and the accessions supporting each call.

Important columns:

- `cluster_id`: VSEARCH cluster id.
- `accession`: sequence accession in the cluster.
- `accession_taxid`: taxid resolved from the accession mapping.
- `assigned_name`: consensus taxon name.
- `assigned_taxid`: consensus taxid.
- `assigned_rank`: rank of the consensus call.
- `disambiguation`: taxa represented inside the cluster.

If a cluster contains accessions from multiple taxa, the barcode sequence is not unique for those taxa at the selected clustering threshold.

## Parsed OBI Table

`parsed_obipcr/` contains one TSV per primer/database pair when OBI is used. It includes:

- binding coordinates
- amplicon length
- primer and match sequences
- mismatch positions from the primer and genome perspective
- mismatch severity
- GC and Tm metrics
- dimer/hairpin indicators
- hit counts and length spread per source sequence

## Divergence Tables

Species-divergence outputs are advisory. The main table contains:

- `accession`
- `species_taxid`
- `subcluster`
- `max_intra_id`
- `mean_intra_id`
- `flag`

Flags are `core`, `divergent`, or `not_screened`.
