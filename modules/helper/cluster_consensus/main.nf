process CLUSTER_CONSENSUS {

    tag "${meta.primer}|${meta.db}"
    label 'process_medium'

    publishDir "${params.outdir}/consensus", mode: 'copy'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cluster_taxonomy)
    path taxdump

    output:
    tuple val(meta), path("*.cluster_consensus.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
    set -euo pipefail

    python3 ${moduleDir}/cluster_consensus.py \\
        --input ${cluster_taxonomy} \\
        --taxdump ${taxdump} \\
        --output ${prefix}.cluster_consensus.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        taxidtools: \$(python3 -c "import taxidTools; print(taxidTools.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
