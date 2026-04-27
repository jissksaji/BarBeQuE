process COMPUTE_BUFFER {

    tag "${meta.id}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    env buffersize
    path "versions.yml", emit: versions

    script:
    """
    buffersize=\$(seqkit stats ${db} | awk 'NR>1 {print \$NF}' | tr -d ',' | sort -nr | head -n 1)
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit //')
    END_VERSIONS
    """
    }