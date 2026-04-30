process VSEARCH_DEREPLICATION {

    tag "${meta.primer}|${meta.db}"

    label 'short_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.derep.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.primer}_${meta.db}"
    """
    vsearch \\
    --fastx_uniques ${fasta} \\
    --fastaout ${prefix}.derep.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1 | sed 's/vsearch v//;s/,.*//')
    END_VERSIONS
    """
}
