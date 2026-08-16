# Usage

## Normal Benchmarking

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --input primers.tsv \
  --dbs refseq_mito,midori_co1 \
  --reference_base /path/to/references \
  --run_name primer_benchmark \
  --outdir results
```

`--run_name` is required when `--dbs` is supplied.

## Primer Inputs

### Samplesheet

Use a tab-separated file:

```text
primer  fwd  rev  min  max
COI     GGWACWGGWTGAACWGTWTAYCCYCC  TAIACYTCIGGRTGICCRAARAAYCA  100  500
```

Column names must be:

- `primer`
- `fwd`
- `rev`
- `min`
- `max`

### Primer FASTA Directory

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --input primer_fastas/ \
  --primer_min 100 \
  --primer_max 500 \
  --dbs refseq_mito \
  --reference_base /path/to/references \
  --run_name fasta_primer_test
```

Each `.fa`, `.fasta`, or `.fna` file is treated as one primer set. Headers with `fwd`, `forward`, `rev`, or `reverse` are collapsed by direction into consensus primers. A plain two-record FASTA without direction labels is treated as forward then reverse.

### Named Primer Sets

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --primer_set COI,Fish16S \
  --dbs midori_co1,midori_lrrna \
  --reference_base /path/to/references \
  --run_name catalog_primers
```

Named sets come from the FooDMe2 primer catalog.

## Database Inputs

### Installed Databases

```bash
--dbs refseq_mito,midori_co1
```

Valid ids are defined in `conf/resources.config`. Use `--list_dbs` to print them.

### Custom Database

```bash
--custom_db /path/to/custom.fasta
```

Optional cleaning:

```bash
--custom_db_filter \
--custom_db_filter_pattern "pattern_to_remove" \
--custom_db_min_length 100 \
--custom_db_max_length 500 \
--custom_db_max_n 0
```

## Taxonomy Options

Use an installed reference base:

```bash
--reference_base /path/to/references
```

Or pass files explicitly:

```bash
--taxdump /path/to/new_taxdump \
--accession_taxonomy /path/to/nucl_gb.accession2taxid
```

Restrict a selected database to a taxon and its descendants before analysis:

```bash
--taxid 9606
```

Remove records assigned to taxids in FooDMe2's built-in blocklist:

```bash
--blocklist
```

This option is off by default. When enabled, BarBeQuE downloads the blocklist
from a pinned FooDMe2 revision, filters each selected FASTA before in-silico
PCR, and writes the filtered database and removal counts to
`blocklist_filtered/`. An accession-to-taxid mapping is required.

## In-Silico PCR Options

Default:

```bash
--insilico_tool obipcr
--obipcr_mismatches 2
--obipcr_fixed_3prime 3
```

`--obipcr_fixed_3prime` forbids mismatches in the last N bases of each primer,
mirroring a real PCR: a mismatch at the 3' end stops the polymerase from
extending, while 5' mismatches are still tolerated up to `--obipcr_mismatches`.
It is applied by appending `#` to those bases in the primer handed to obipcr,
e.g. `GGGCAATCCTGAGCCAA` becomes `GGGCAATCCTGAGCC#A#A#`. Set it to `0` to allow
mismatches anywhere in the primer.

Alternative:

```bash
--insilico_tool cutadapt
--cutadapt_mismatches 2
```

## Clustering Options

```bash
--cluster_id 0.97
```

This controls the identity threshold used by `vsearch --cluster_fast`.

## Masking

```bash
--mask --read_length 150
```

For single-end data:

```bash
--mask --single_end --read_length 400
```

## Divergence Screening And Exclusion

Advisory screening:

```bash
--screen_species_divergence \
--species_divergence_id 0.99 \
--species_hclust_method average \
--species_hclust_max_seqs 2000
```

Curated filtering before clustering:

```bash
--exclude_accessions bad_accessions.txt
```

`bad_accessions.txt` contains one accession per line. Blank lines and lines starting with `#` are ignored. Matching is version-insensitive.

## Target-Taxon Reports

```bash
--taxon Mammalia
```

This adds taxon-focused coverage outputs.

```bash
--completeness_table
```

This adds database completeness summaries using the value supplied to `--taxon`.

## Interactive Dashboard

```bash
--interactive
```

This launches the Streamlit dashboard after the analysis finishes.

## Standalone Hierarchical Clustering

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --hierarchical_clustering \
  --custom_db custom.fasta \
  --reference_base /path/to/references \
  --run_name custom_screen \
  --outdir results_hclust
```

This mode screens a custom FASTA directly. Do not provide `--input`, `--primer_set`, or `--dbs`.
