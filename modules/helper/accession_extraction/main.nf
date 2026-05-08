process ACCESSION_EXTRACTION {

    tag "${db.meta}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    path ("*.accessions.txt"), emit: txt
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${db}"

    """
    awk '/^>/ {
        sub(/^>/, "", \$1)
        split(\$1, a, ".")
        print a[1]
    }' ${db} \\
    | sort -u \\
    > ${prefix}.accessions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """
}
