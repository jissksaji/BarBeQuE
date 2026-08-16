process BLASTN_CORE_NT {

    tag "${meta.id}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1' :
        'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta), path(query_fasta)
    val(core_nt_db)
    val(taxid)

    output:
    tuple val(meta), path("${meta.id}.core_nt.blastn.tsv"), emit: tsv
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def taxid_arg = taxid ? "-taxids ${taxid}" : ''
    """
    blastn \\
        -query ${query_fasta} \\
        -db ${core_nt_db} \\
        ${taxid_arg} \\
        -num_threads ${task.cpus} \\
        -out ${meta.id}.core_nt.blastn.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
END_VERSIONS
    """
}
