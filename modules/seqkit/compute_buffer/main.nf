process COMPUTE_BUFFER {

    tag "${meta.id}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    env buffersize, emit: buffersize
    path "versions.yml", emit: versions

    script:
    """
    buffersize=\$(seqkit stats -T ${db} | awk -F'\t' 'NR==2 {print \$8 * 5 + 50000}')
    buffersize=\$((buffersize * 100))


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit //')
    END_VERSIONS
    """
}
