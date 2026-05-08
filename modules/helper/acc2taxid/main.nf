process ACC2TAXID {

    tag "${meta.db}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    path ("*.acessions.txt"), emit: txt
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.db}"
}
