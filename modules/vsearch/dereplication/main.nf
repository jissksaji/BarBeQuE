process VSEARCH_DEREPLICATION {
    tag "${meta.primer}|${meta.db}"

    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    input:
    tuple val(meta),path(fasta)

    output: 
    tuple val(meta),path("*.derep.fasta"),   emit: fasta
    path "versions.yml" ,                    emit : versions

    script:
    def prefix = "${meta.primer}_${meta.db}"
    """
    vsearch \\
    --fastx_uniques ${fasta} \\
    --fastaout ${prefix}.derep.fasta \\
    --threads ${task.cpu}
    """
}


