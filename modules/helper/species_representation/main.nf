process SPECIES_REPRESENTATION {

    tag "${meta.primer}|${meta.db}"
    label 'process_low'



    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(coverage_tsv), path(consensus_tsv)

    output:
    tuple val(meta), path("*.representation.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
    species_representation.py \\
        --coverage ${coverage_tsv} \\
        --consensus ${consensus_tsv} \\
        --output ${coverage_tsv.name.replaceAll('.tax_coverage', '')}.representation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
