process PARSE_PRIMERS {

    tag "${meta.id}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(primer_input)

    output:
    tuple val(meta), path("${meta.id}.tsv")     , emit: samplesheet
    path("${meta.id}.primer_warnings.txt")      , emit: warnings
    path("versions.yml")                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_primers.py \\
        --input ${primer_input} \\
        --min ${meta.min} \\
        --max ${meta.max} \\
        --out ${meta.id}.tsv \\
        --warnings ${meta.id}.primer_warnings.txt

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
        python: \$(python3 --version | sed 's/Python //')
END_VERSIONS
    """
}
