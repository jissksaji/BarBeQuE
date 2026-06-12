process COMPLETENESS_TABLE {

    tag "${meta.id}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/completeness", mode: 'copy'

    input:
    tuple val(meta), path(db_taxids)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    completeness_table.py \\
        --species "${params.taxon}" \\
        --taxids_counts ${db_taxids} \\
        --taxdump ${params.taxdump ?: params.references.taxdump} \\
        --output ${prefix}.completeness.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
