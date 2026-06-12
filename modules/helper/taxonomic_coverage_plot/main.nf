process TAXONOMIC_COVERAGE_PLOT {
    debug true

    maxForks 1
    debug true
    cache false

    label 'process_low'
    publishDir "${params.outdir}/taxonomic_coverage/", mode: 'copy'

    conda "${moduleDir}/environment.yml"

    input:
    path tsv_files

    output:
    path '*.html', emit: plot
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: "combined_taxonomic_coverage"
    """
    taxonomic_coverage_plot.py \\
        --input ${tsv_files} \\
        --output "${prefix}.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
