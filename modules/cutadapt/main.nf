process CUTADAPT_INSILICOPCR {

    tag "${meta.primer}|${meta.db}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/cutadapt:5.2--py313h8c92656_1"

    input:
    tuple val(meta), path(db), env(buffersize)

    output:
    tuple val(meta), path('*_insilico.fasta'), emit: fasta
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    def rev_r = meta.rev.reverse()

    """
cutadapt ${args} \
    -g "${meta.fwd}...${rev_r}" \
    --overlap ${params.cutadapt_overlap} \
    -e ${params.cutadapt_error_rate} \
    --cores ${task.cpus} \
    --buffer-size \$buffersize \
    --discard-untrimmed \
    --minimum-length ${meta.min} \
    --maximum-length ${meta.max} \
    -o "${prefix}_insilico.fasta" \
    --no-indels \
    --revcomp \
    ${db}

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
        cutadapt: \$(cutadapt --version)
END_VERSIONS
    """
}
