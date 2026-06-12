process DB_DISTRIBUTION {

    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/db_distribution", mode: 'copy'

    input:
    tuple val(meta), path(taxids_counts)

    output:
    tuple val(meta), path("*.db_distribution.tsv"), emit: distribution
    tuple val(meta), path("*.db_distribution.html"), emit: plot
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    db_distribution.py \
        --input ${taxids_counts} \
        --output ${prefix}.db_distribution.tsv

    db_distribution_plot.py \
        --input  ${prefix}.db_distribution.tsv \
        --output ${prefix}.db_distribution.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        ete3: \$(python -c "import ete3; print(ete3.__version__)")
    END_VERSIONS
    """
}
