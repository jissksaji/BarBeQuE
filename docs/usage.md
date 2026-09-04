# Usage

## Before You Start

Two things have to be settled before the first run.

### 1. Choose exactly one primer input

`--input` and `--primer_set` are mutually exclusive, and one of them is required.
Passing both, or neither, stops the run straight away.

```bash
--input primers.tsv        # a samplesheet, a primer FASTA, or a directory of them
--primer_set COI,Fish16S   # named sets from the FooDMe2 catalog
```

### 2. Make the reference data available

Either install everything once with `--build_references`, which lays out the
databases, the NCBI taxdump and the accession-to-taxid mapping under a single
directory:

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --build_references \
  --reference_base /path/to/references
```

Every later run then needs only `--reference_base /path/to/references`. See
[Reference Installation](build_references.md).

Or skip the reference base and point at each piece yourself, which is the usual
route when you are benchmarking your own database:

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --input primers.tsv \
  --custom_db /path/to/custom.fasta \
  --taxdump /path/to/new_taxdump \
  --accession_taxonomy /path/to/nucl_gb.accession2taxid \
  --outdir results
```

The three go together: `--custom_db` supplies the sequences, `--taxdump` the NCBI
taxonomy tree, and `--accession_taxonomy` the accession-to-taxid mapping that links
the two. Without `--reference_base` there is nothing for the pipeline to fall back
on, so a missing one is an error rather than a default.

Note that `--custom_db` takes precedence: if you pass both `--custom_db` and
`--dbs`, only the custom database is used.

## Seeing What Is Available

Both of these print a list and exit without running an analysis, so they are safe
to try first.

Installed databases:

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --list_dbs
```

Each line gives the id you pass to `--dbs` and where the data came from.

Named primer sets:

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --list_primers
```

This prints the catalog as JSON: the name to pass to `--primer_set`, a description,
the target markers, the sequencing platform, and the publication DOI. It needs
internet access.

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

`--input` accepts a single FASTA file or a directory of them. A single file is
recognised by its **content** - a first non-blank line starting with `>` - so its
extension does not matter and a FASTA named `primers.txt` still works. Inside a
directory, only `.fa`, `.fasta` and `.fna` files are read and everything else is
ignored.

Both forms require global amplicon bounds, since a FASTA carries no length
information:

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

Named sets are pulled at run time from the
[FooDMe2](https://github.com/bio-raum/FooDMe2) primer catalog, so any primer system
published there can be benchmarked here without writing a samplesheet of your own.
Run `--list_primers` to see what is on offer, then pass one or more names as a
comma-separated list.

The catalog and its primer FASTAs are fetched over HTTPS from a FooDMe2 revision
pinned inside BarBeQuE, so a given pipeline version always tests the same sequences
even if FooDMe2 changes upstream. This step needs internet access.

Each set brings its own amplicon `min`/`max` from the catalog, so `--primer_min` and
`--primer_max` are not needed. Downloaded sets go through the same parsing step as a
primer FASTA and end up as the same `primer / fwd / rev / min / max` rows as a
samplesheet, so nothing downstream behaves differently.

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

The same optional cleaning can be applied to custom or installed databases:

```bash
--db_filter \
--db_filter_pattern "pattern_to_remove" \
--db_filter_min_length 100 \
--db_filter_max_length 500 \
--db_filter_max_n 0
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
--mask
```

Masking rewrites each in-silico amplicon to look like what a sequencer would
actually deliver, so a barcode is only credited with the part real reads would
cover. Without it, every amplicon is treated as if it were sequenced end to end.

**Paired-end is the default.** Each amplicon keeps `--read_length` bases at both
ends and has its middle replaced by 20 `N`s. Amplicons short enough for the two
reads to overlap - length up to `2 x read_length` - are left untouched.

```bash
--mask --read_length 150
```

**Single-end** truncates each amplicon to `--read_length` bases and inserts no
`N`s, since there is no second read coming back from the far end:

```bash
--mask --single_end --read_length 400
```

`--read_length` defaults to **150 for paired-end** and **400 for single-end**,
matching common 2x150 and single-400 chemistries. Both examples above therefore
just restate the default - pass `--read_length` only when modelling a different
read length.

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

## Interactive Dashboard

```bash
--interactive
```

Once the analysis has finished, this starts a Streamlit dashboard on port **8501**
of the machine running the pipeline, reachable at <http://localhost:8501>. Any
Streamlit instance already serving this app is stopped first, so the port has to be
free for BarBeQuE to use it.

The dashboard is **experimental**. The main page is the most dependable one; several
of the pages below it were written around particular analyses and may error out or
come up empty for your run. Treat it as an extra view of the results rather than the
results themselves - everything it draws is read from files already written to
`--outdir`.

## Next Steps

Once a run finishes, see [Outputs](output.md) for what lands in `--outdir` and how to
read it. If a run fails, [Troubleshooting](troubleshooting.md) covers the common
causes.
