process DB_FILTER {

    tag "${meta.id}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.cleaned.fasta"), emit: fasta
    path("versions.yml"),                       emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pattern = params.db_filter_pattern

    // Length and ambiguity filtering is enabled when both length bounds are set.
    def do_length_filter = params.db_filter_min_length != null && params.db_filter_max_length != null
    def length_filter = do_length_filter
        ? "| seqkit seq" +
          " --min-len ${params.db_filter_min_length}" +
          " --max-len ${params.db_filter_max_length}" +
          " --threads ${task.cpus}"
        : ""
    // SeqKit has no --max-ambig option. Reject records containing more than
    // the allowed number of IUPAC ambiguity symbols with a sequence regex.
    def ambiguity_filter = do_length_filter && params.db_filter_max_n != null
        ? "| seqkit grep --by-seq --use-regexp --ignore-case --only-positive-strand --invert-match" +
          " --pattern \"[RYKMSWBDHVN]([^RYKMSWBDHVN]*[RYKMSWBDHVN]){${params.db_filter_max_n}}\"" +
          " --threads ${task.cpus}"
        : ""

    """
    seqkit grep -n -v -r -i \
        --threads ${task.cpus} \
        -p "${pattern}" \
        ${db} \
        ${length_filter} \
        ${ambiguity_filter} \
        -o ${prefix}.cleaned.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit //')
    END_VERSIONS
    """
}
