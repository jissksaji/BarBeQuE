process CUTADAPT_INSILICOPCR {

    tag "${meta.primer}|${meta.db}"

    label 'medium_parallel'

    // conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/cutadapt:5.2--py313h8c92656_1"

    input:
    tuple val(meta), path(db), env('buffersize')

    output:
    tuple val(meta), path('*_insilico.fasta'), emit: fasta
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    def rev_r = meta.rev.toUpperCase().reverse().tr('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')

    def mismatches = params.cutadapt_mismatches as int

    def fwd_len = meta.fwd.length()
    def rev_len = rev_r.length()

    // mathematically exact overlap formula to guarantee 'mismatches' allowed errors
    def e_value = mismatches == 0 ? 0 : mismatches + 0.5
    def fwd_overlap = mismatches == 0 ? fwd_len : Math.ceil((mismatches * fwd_len) / e_value) as int
    def rev_overlap = mismatches == 0 ? rev_len : Math.ceil((mismatches * rev_len) / e_value) as int

    """
cutadapt ${args} \
    -g "${meta.fwd};e=${e_value};o=${fwd_overlap}...${rev_r};e=${e_value};o=${rev_overlap}" \
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
