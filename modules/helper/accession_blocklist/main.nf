process ACCESSION_BLOCKLIST {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta), path(tsv)
    path accession_blocklist

    output:
    tuple val(meta), path("*.accession_filtered.fasta"), emit: fasta
    tuple val(meta), path("*.accession_filtered.tsv"), emit: tsv
    tuple val(meta), path("*.accession_blocklist_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
    filter_accession_blocklist.py \
        --fasta ${fasta} \
        --tsv ${tsv} \
        --accession-blocklist ${accession_blocklist} \
        --fasta-output ${prefix}.accession_filtered.fasta \
        --tsv-output ${prefix}.accession_filtered.tsv \
        --summary ${prefix}.accession_blocklist_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
