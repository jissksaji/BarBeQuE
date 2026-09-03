# Outputs

Outputs are written under `--outdir`, default `results`.

## Core Directories

| Directory | Contents |
| --- | --- |
| `pipeline_info/` | Staged samplesheet and software version metadata. |
| `reports/` | MultiQC reports, one per primer/database group. |
| `raw/obipcr/` | Raw OBI in-silico PCR output. |
| `parsed_obipcr/` | Parsed OBI results with primer mismatch and amplicon metrics. |
| `primers/` | The samplesheet generated from primer FASTA input, plus any parser warnings. |
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
| `accession_blocklist/` | `--accession_blocklist` | Per-primer/database summaries of listed, matched, and removed accessions. |
| `tax_coverage/` | `--taxon` | Taxon-focused coverage and species representation tables. |

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
Consensus assignments use `--consensus_fraction` (default `1.0`). The fraction is
calculated over accessions with valid taxonomy; accessions without a usable taxid
do not vote and are not included in the denominator.

## Parsed OBI Table

`parsed_obipcr/` contains one TSV per primer/database pair when OBI is used. When
`--accession_blocklist` is supplied, these are the filtered tables. They include:

- binding coordinates
- amplicon length
- primer and match sequences
- mismatch positions from the primer and genome perspective
- mismatch severity
- GC and Tm metrics
- dimer/hairpin indicators
- hit counts and length spread per source sequence

## Accession Blocklist Summary

`accession_blocklist/` contains one `*.accession_blocklist_summary.tsv` per
primer/database pair:

- `listed_accessions`, `matched_accessions`, and `unmatched_accessions`
- input, kept, and removed FASTA record counts
- input, kept, and removed parsed TSV row counts
