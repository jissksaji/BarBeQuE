process COVERAGE_ANALYSIS {

    tag "${meta.primer}|${name}"
    label 'short_parallel'

    publishDir "${params.outdir}/coverage_analysis/blast", mode: 'copy', pattern: '*.blast.tsv'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(centroids), val(name), path(fasta), path(db)

    output:
    path "*.coverage.tsv", emit: tsv
    path "*.blast.tsv", emit: blast
    path "versions.yml", emit: versions

    script:
    // Split the file name into its parts for the report
    // Example: 155109_Cyclopia_subternata_species -> 155109, Cyclopia subternata, species
    def parts = name.tokenize('_')
    def taxid = parts.first()
    def rank = parts.last()
    def taxon_name = parts[1..-2].join(' ')
    """
    # Record the tool versions first, because the script can stop early below

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n 1 | sed 's/blastn: //')
    END_VERSIONS

    # Start the report with a header line

    echo -e "primer\\ttaxid\\tname\\trank\\trecords\\twith_region" > ${meta.primer}_${name}.coverage.tsv

    # Start the BLAST table with a header line (it stays like this for an empty taxid file)

    echo -e "centroid\\trecord\\tidentity\\tcentroid_cov\\taligned\\tcentroid_len\\trecord_from\\trecord_to" > ${meta.primer}_${name}.blast.tsv

    # Count the records in the taxid file

    records=\$(grep -c '>' ${fasta} || true)

    # Empty taxid file: nothing to search, so the report row is 0 and 0

    if [ "\$records" -eq 0 ]; then
        echo -e "${meta.primer}\\t${taxid}\\t${taxon_name}\\t${rank}\\t0\\t0" >> ${meta.primer}_${name}.coverage.tsv
        exit 0
    fi

    # Search the centroids against the taxid records
    # A record has the region if one hit covers most of a centroid at high identity
    # -max_target_seqs: BLAST stops at 500 hits per centroid by default, so allow many more
    # -dust no: do not mask simple repeats like ATATAT
    # Each line: which centroid hit which record, how well, and where in the record

    blastn -query ${centroids} -db db/${name} \\
        -perc_identity ${params.coverage_identity} \\
        -qcov_hsp_perc ${params.coverage_query_cov} \\
        -max_target_seqs 1000000 \\
        -dust no \\
        -num_threads ${task.cpus} \\
        -outfmt "6 qseqid sseqid pident qcovhsp length qlen sstart send" >> ${meta.primer}_${name}.blast.tsv

    # Count each record once, even if several centroids hit it
    # tail skips the header line, cut -f2 keeps the record column
    
    with_region=\$(tail -n +2 ${meta.primer}_${name}.blast.tsv | cut -f2 | sort -u | wc -l)

    echo -e "${meta.primer}\\t${taxid}\\t${taxon_name}\\t${rank}\\t\$records\\t\$with_region" >> ${meta.primer}_${name}.coverage.tsv
    """
}
