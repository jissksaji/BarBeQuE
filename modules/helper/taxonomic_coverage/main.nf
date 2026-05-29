process TAXONOMIC_COVERAGE {

    maxForks 1
    debug true
    cache false

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(clusters), path(db_taxids)
    val taxonomy

    output:
    tuple val(meta), path('*.tsv'), emit: tsv
    tuple val(meta), path('*.nwk'), emit: nwk
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}_--${taxonomy}--"
    """

    ete.py --taxon "${taxonomy}" \\
    --reference ${db_taxids} \\
    --report ${clusters} \\
    --output ${prefix}.tax_coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
