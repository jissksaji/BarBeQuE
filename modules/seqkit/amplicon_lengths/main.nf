process AMPLICON_LENGTH {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/amplicon_lengths", mode: 'copy'


    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml",            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
    seqkit fx2tab -l -n ${fasta} > ${prefix}.amplicon_lengths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}