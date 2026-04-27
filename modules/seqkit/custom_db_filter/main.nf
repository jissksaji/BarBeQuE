process CUSTOM_DB_FILTER {

    tag "${meta.id}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.cleaned.fasta"), emit: fasta
    //path("*.stats.txt"),                      emit: stats
    path("versions.yml"),                     emit: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    PATTERN="environmental|uncultured|metagenom|unidentified|unknown|unclassified|incertae[ ._-]?sedis|synthetic|artificial|construct|recombinant|transgenic|predicted|XM_|XR_|unverified|low[ ._-]?quality|contaminant|chimeric|misidentified"

    seqkit grep -n -v -r -i \\
        --threads ${task.cpus} \\
        -p "\$PATTERN" \\
        ${db} \\
        -o ${prefix}.cleaned.fasta

    #seqkit stats ${db} ${prefix}.cleaned.fasta \\
    #    --tabular \\
    #    -o ${prefix}.stats.txt 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit //')
    END_VERSIONS
    """
}