process BLOCKLIST_FILTER {

    tag "${meta.id}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)
    path accession_taxonomy
    path blocklist

    output:
    tuple val(meta), path("*.blocklist_filtered.fasta"), emit: fasta
    tuple val(meta), path("*.blocklist_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_blocklist.py \\
        --fasta ${db} \\
        --accession-taxid ${accession_taxonomy} \\
        --blocklist ${blocklist} \\
        --output ${prefix}.blocklist_filtered.fasta \\
        --summary ${prefix}.blocklist_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
