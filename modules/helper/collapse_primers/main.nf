process COLLAPSE_PRIMERS {

    tag "${meta.id}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "python:3.14-slim"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.collapsed.fasta"), emit: fasta
    path("versions.yml"), emit: versions

    script:
    """
    collapse_primers.py --fasta ${fasta} --prefix ${meta.id} --out ${meta.id}.collapsed.fasta

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
        python: \$(python3 --version | sed 's/Python //')
END_VERSIONS
    """
}
