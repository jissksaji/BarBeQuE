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

### Primer FASTA

`--input` accepts a single `.fa`/`.fasta`/`.fna` file or a directory of them. Both require global
amplicon bounds, since a FASTA carries no length information:

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

Each file becomes one or more primer sets, named after the filename. Records are grouped by the
prefix of their header, same-length variants of a primer are collapsed into one IUPAC-degenerate
sequence, and everything is validated before the run starts. The generated samplesheet is published
to `primers/` so you can check exactly what was benchmarked.

See [Primer Input](primer_input.md) for the full rules: prefix grouping, collapsing and splitting,
set naming, and what makes a file invalid.

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

## In-Silico PCR Options

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

## Clustering Options

```bash
--cluster_id 0.97
--consensus_fraction 0.8
```

`--cluster_id` controls the identity threshold passed to `vsearch --cluster_fast`
as `--id`. It accepts a value from `0.0` to `1.0` and defaults to `0.97`; provide
another value on the Nextflow command line to override it for a run.
`--consensus_fraction` controls the minimum fraction of taxonomically assigned
accessions in a cluster that must support a taxon. It must be greater than `0.5`
and at most `1.0`. The default, `1.0`, is a strict lowest-common-ancestor call;
lower values allow a supported majority assignment despite minority outliers.

## Masking

```bash
--mask --read_length 150
```

For single-end data:

```bash
--mask --single_end --read_length 400
```

## Excluding Unwanted Accessions

```bash
--accession_blocklist unwanted_accessions.txt
```

The file holds one accession per line. Blank lines and text after `#` are
ignored. Matching is case- and version-insensitive, so either a base or
versioned accession can be supplied:

```text
# misannotated records
MK123456
AY846379.1
```

Filtering happens immediately after OBI-PCR parsing. Matching rows are removed
from the parsed TSV and the amplicon FASTA, so they cannot reach masking,
amplicon summaries, clustering, taxonomy, consensus, or reports. Audit counts
are written to `accession_blocklist/*.accession_blocklist_summary.tsv`.

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
