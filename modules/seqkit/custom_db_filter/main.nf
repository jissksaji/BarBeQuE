process CUSTOM_DB_FILTER {

    //removes the sequences like environmental or  unclassified or unknown  from the db
    //works purley by grepping and removing the sequnces which includes these terms from the header
    //nextflow.config has the defaults

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
    def pattern = params.custom_db_filter_pattern
    
    //length filter if both min and max are given
    def do_length_filter =params.custom_db_min_length && params.custom_db_max_length
    def len_filter = do_length_filter ?
         "|seqkit seq" + 
         " --min-len ${params.custom_db_min_length}"+
         " --max-len ${params.custom_db_max_length}"+
        (params.custom_db_max_n ? " --max-ambig ${params.custom_db_max_n}": "")+
        " --threads ${task.cpus}" : ""


    """

    seqkit grep -n -v -r -i \\
        --threads ${task.cpus} \\
        -p "${pattern}" \\
        ${db} \\
        ${len_filter} \\
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