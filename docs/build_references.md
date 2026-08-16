# Reference Installation

`--build_references` is an installation workflow. It does not run primer benchmarking.

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --build_references \
  --reference_base /path/to/references
```

References are installed under:

```text
<reference_base>/barbeque/<reference_version>/
```

`--reference_version` defaults to `1.1`.

## Installed Assets

The workflow installs three kinds of data.

### Reference FASTAs

All installable databases are declared in `conf/resources.config`.

- `refseq_mito`
- `refseq_plastid`
- `refseq_plasmid`
- `midori_lrrna`
- `midori_srrna`
- `midori_cytb`
- `midori_co1`
- `midori_co2`
- `midori_co3`
- `mitofish`
- `metafish`
- `silva_ssu`

`metafish` is downloaded from the MetaFish library CSV and converted to FASTA during installation.

`silva_ssu` contains the combined SILVA SSU collection (16S and 18S). It is
installed once instead of publishing the same source under two database IDs.
The workflow normalizes SILVA-style accession headers with extra numeric
suffixes during taxonomy lookup.

`core_nt` remains a special case. It is a preformatted BLAST database and is not used as an in-silico-PCR FASTA through `--dbs`.

Every installed sequence database has a `.fasta` filename. To add another
database, copy an existing block in `conf/resources.config`:

```groovy
'example' {
  urls = ['https://example.org/release/example.fasta.gz']
  format = 'fasta'
  release = '2026-01-15'
  db = "${params.reference_base}/barbeque/${params.reference_version}/databases/example/example.fasta"
  description = 'Example database'
}
```

The build workflow loops over these entries and installs all of them. Each
download task records the release and final SHA-256 checksum.

### Primer FASTAs

Primer FASTAs are resolved from the FooDMe2 catalog by
`lib/PrimerCatalog.groovy`. The FooDMe2 Git revision is pinned in that class so
the catalog and its FASTAs remain reproducible. A catalog-wide failure stops the
installation.

### Taxonomy Files

Taxonomy support is enabled by default. `--install_taxdump` installs both:

- NCBI `new_taxdump`
- NCBI `nucl_gb.accession2taxid`

Both are needed by most analysis runs unless equivalent files are supplied later with `--taxdump` and `--accession_taxonomy`.

They are installed together under:

```text
<reference_base>/barbeque/<reference_version>/taxonomy/
```

## MIDORI Version

`--midori_version` controls which MIDORI release is used for configured MIDORI downloads. Keep this pinned for reproducible analyses.

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --build_references \
  --reference_base /path/to/references \
  --midori_version 271_2026-04-07
```

## Completion Behavior

On success, all final references are available below `reference_base`.
