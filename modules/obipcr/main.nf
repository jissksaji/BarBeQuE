process OBIPCR_INSILICOPCR {

    tag "${meta.primer}|${meta.db}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path('*_insilico.fasta'), emit: fasta
    tuple val(meta), path('*_raw.fasta'), emit: raw_fasta, optional: true
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    def mismatches = params.containsKey('obipcr_mismatches')
        ? params.obipcr_mismatches as int
        : params.cutadapt_mismatches as int

    def fwd_primer = meta.fwd
    def rev_primer = meta.rev

    """
    obipcr ${args} \\
        --forward "${fwd_primer}" \\
        --reverse "${rev_primer}" \\
        --allowed-mismatches ${mismatches} \\
        --min-length ${meta.min} \\
        --max-length ${meta.max} \\
        --max-cpu ${task.cpus} \\
        --fasta-output \\
        --no-progressbar \\
        --skip-empty \\
        --only-complete-flanking \\
        "${db}" \\
        | tee "${prefix}_raw.fasta" \\
        | seqkit replace -p '_sub\\[.*' -r '' \\
        | seqkit seq -w 0 -u \\
        > "${prefix}_insilico.fasta"

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    obipcr: "\$(obipcr --version 2>&1 | head -n 1)"
END_VERSIONS
    """
}
