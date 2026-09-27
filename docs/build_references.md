# Reference Installation

`--build_references` installs the reference data and exits. It doesn't run an analysis.

```bash
nextflow run main.nf \
  -profile conda \
  --build_references \
  --reference_base /path/to/references
```

Everything is installed under `<reference_base>/barbeque/<reference_version>/` (the default
`--reference_version` is `1.1`):

| Folder | Contents |
| --- | --- |
| `databases/<id>/` | Every database in `conf/resources.config` |
| `primers/` | Primer FASTAs from the FooDMe2 catalog |
| `taxonomy/` | NCBI `new_taxdump` and `nucl_gb.accession2taxid` |

## Databases

`refseq_mito`, `refseq_plastid`, `refseq_plasmid`, `midori_lrrna`, `midori_srrna`, `midori_cytb`,
`midori_co1`, `midori_co2`, `midori_co3`, `mitofish`, `metafish`, `silva_ssu`, `silva_lsu`,
`its2_global`

`--list_dbs` prints the same list. `core_nt` is also listed there, but it isn't downloaded:
it's an externally managed BLAST database.

`--midori_version` (default `271_2026-04-07`) sets which MIDORI release is downloaded.

To add a database, copy an existing entry in `conf/resources.config`:

```groovy
'example' {
  urls = ['https://example.org/example.fasta.gz']
  format = 'fasta'
  release = '2026-01-15'
  db = "${params.reference_base}/barbeque/${params.reference_version}/databases/example/example.fasta"
  description = 'Example database'
}
```

## Taxonomy

`--install_taxdump` (on by default) installs the NCBI taxdump and the accession-to-taxid mapping.
Most runs need both. Set it to `false` if you'll supply your own later with `--taxdump` and
`--accession_taxonomy`.
