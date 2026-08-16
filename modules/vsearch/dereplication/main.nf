process VSEARCH_DEREPLICATION {

    tag "${meta.primer}|${meta.db}|sp${meta.species}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.derep.fasta"), emit: fasta
    tuple val(meta), path("*.derep.uc"), emit: uc
    path "versions.yml", emit: versions

    script:
    def prefix = meta.species ? "${meta.primer}_${meta.db}.sp_${meta.species}" : "${meta.primer}_${meta.db}"
    """
    vsearch \\
    --fastx_uniques ${fasta} \\
    --sizeout \\
    --uc ${prefix}.derep.uc \\
    --fastaout ${prefix}.derep.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n1 | sed 's/vsearch v//;s/,.*//')
    END_VERSIONS
    """
}
