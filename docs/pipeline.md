# Pipeline Workflow

`main.nf` is the routing layer for BarBeQuE. It enables DSL2, prints parameter help through `nf-schema`, validates parameters through `WorkflowMain` and `WorkflowPipeline`, and then chooses one of three workflows.

## Entry Points

```groovy
if (params.build_references) {
    BUILD_REFERENCES()
} else if (params.hierarchical_clustering) {
    HIERARCHICAL_CLUSTERING_WORKFLOW()
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
3. Optional `--blocklist` removes records assigned to FooDMe2-blocklisted taxids.
4. Optional `--custom_db_filter` cleans a custom database.
5. `BARBEQUE` resolves primers.
6. In-silico PCR runs with `obipcr` or `cutadapt`.
7. Amplicons can be masked, length-profiled, taxonomically mapped, and clustered.
8. Cluster membership is joined to accession taxonomy.
9. Consensus taxonomy, database distribution, optional taxon coverage, optional completeness, and MultiQC reports are written.

See [barbeque.md](barbeque.md) for the step-by-step analysis workflow.

## Reference Installation Mode

`--build_references` skips analysis and installs reference assets under:

```text
<reference_base>/barbeque/<reference_version>/
```

It downloads configured FASTA databases, primer FASTAs from the FooDMe2 catalog, and optionally NCBI taxonomy/accession mapping files. See [build_references.md](build_references.md).

## Standalone Hierarchical Clustering Mode

`--hierarchical_clustering` uses `--custom_db` directly as the screened sequence set. It does not run primer benchmarking.

It requires:

- `--custom_db`
- either `--reference_base`, or both `--taxdump` and `--accession_taxonomy`

It rejects `--input`, `--primer_set`, and `--dbs` because primers and selected databases are not part of this mode.

## Completion

`PIPELINE_COMPLETION` always runs after the selected workflow. It handles final pipeline bookkeeping and report artefacts.
