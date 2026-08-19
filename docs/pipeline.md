# Pipeline Workflow

`main.nf` is the routing layer for BarBeQuE. It enables DSL2, prints parameter help through `nf-schema`, validates parameters through `WorkflowMain` and `WorkflowPipeline`, and then chooses one of two workflows.

## Entry Points

```groovy
if (params.build_references) {
    BUILD_REFERENCES()
} else {
    DATABASE()
    BARBEQUE(DATABASE.out.db, DATABASE.out.versions)
    if (params.interactive) {
        BARBEQUE.out.consensus.collect() | map { "${params.outdir}" } | INTERACTIVE_RESULTS
    }
}
PIPELINE_COMPLETION()
```

## Normal Benchmarking Mode

This is the default mode. It requires either `--input` or `--primer_set`, plus either `--dbs` or `--custom_db`.

High-level flow:

1. `DATABASE` resolves selected reference FASTAs.
2. Optional `--taxid` filters each database to one taxon and its descendants.
3. Optional `--custom_db_filter` cleans a custom database.
4. `BARBEQUE` resolves primers. `--primer_set` downloads and `--input` FASTAs both pass through
   `PARSE_PRIMERS` into a samplesheet, which `INPUT_CHECK` then validates - see
   [primer_input.md](primer_input.md).
5. In-silico PCR runs with `obipcr` and its output is parsed.
6. Optional `--accession_blocklist` removes matching parsed hits and amplicons.
7. Retained amplicons can be masked, length-profiled, taxonomically mapped, and clustered.
8. Cluster membership is joined to accession taxonomy.
9. Consensus taxonomy, database distribution, optional taxon coverage, optional completeness, and MultiQC reports are written.

See [barbeque.md](barbeque.md) for the step-by-step analysis workflow.

## Reference Installation Mode

`--build_references` skips analysis and installs reference assets under:

```text
<reference_base>/barbeque/<reference_version>/
```

It downloads configured FASTA databases, primer FASTAs from the FooDMe2 catalog, and optionally NCBI taxonomy/accession mapping files. See [build_references.md](build_references.md).

## Completion

`PIPELINE_COMPLETION` always runs after the selected workflow. It handles final pipeline bookkeeping and report artefacts.
