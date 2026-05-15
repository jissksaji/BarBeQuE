process PARSE_UC {

    tag "${meta.primer}|${meta.db}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(uc)

    output:
    tuple val(meta), path("*.cluster_accessions.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
    set -euo pipefail

    awk -F'\\t' '
        BEGIN { OFS = "\\t" }
        \$1 == "S" || \$1 == "H" {
            acc = \$9

            #  cleanup of accession field
            sub(/;size=[0-9]+;?\$/, "", acc)     # strip vsearch size annotation
            sub(/\\.[0-9]+\$/, "", acc)          # strip version (.1, .2, etc.)

            #skip empty or *
            if (acc != "*" && acc != "") {
                print \$2, acc
            }
        }
    ' "${uc}" \\
    | sort -k1,1n -k2,2 \\
    > "${prefix}.cluster_accessions.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
