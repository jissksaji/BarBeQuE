process PARSE_UC {

    tag "${meta.primer}|${meta.db}"
    label 'short_serial'

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
            # Two rules for this awk block, both learned the hard way:
            #   1. escape every dollar sign as \\\$ -- an unescaped one (above all a
            #      dollar immediately followed by a slash) is eaten by Nextflow
            #      interpolation and silently corrupts the generated script
            #   2. no apostrophes, not even in comments -- they close the shell quote
            sub(/;.*/, "", acc)                  # drop everything from the first semicolon:
                                                 # covers SINTAX headers (ACC;tax=k:...,s:...;)
                                                 # and the vsearch size annotation (;size=N;)
            sub(/(\\.[0-9]+)+\$/, "", acc)       # strip version/range (.1, .1.1791, etc.)

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
