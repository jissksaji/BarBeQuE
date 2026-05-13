process TAXONKIT_LCA {

    tag "${meta.primer}|${meta.db}"

    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cluster_taxids)
    path taxdump

    output:
    tuple val(meta), path("*.cluster_lca.tsv"), emit: lca
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    set -euo pipefail

    taxonkit lca \\
        --data-dir "${taxdump}" \\
        -i 2 \\
        -s " " \\
        -D \\
        -U \\
        "${cluster_taxids}" \\
    > "${prefix}.cluster_lca.raw.tsv"

    {
        echo -e "cluster_id\\ttaxids_in_cluster\\tlca_taxid\\tlca_lineage\\tlca_name\\tlca_rank"

        taxonkit lineage \\
            --data-dir "${taxdump}" \\
            -i 3 \\
            -n \\
            -r \\
            "${prefix}.cluster_lca.raw.tsv"
    } > "${prefix}.cluster_lca.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$(taxonkit version 2>&1 | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
