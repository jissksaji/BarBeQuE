process UC_TO_CLUSTER_ACCESSIONS {

    tag "${meta.primer}|${meta.db}"
    label 'medium_parallel'

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

    awk '
        BEGIN { OFS = "\\t" }

        \$1 == "S" || \$1 == "H" {
            cluster = \$2
            acc = \$9
            sub(/^>/, "", acc)
            sub(/\\.[0-9]+\$/, "", acc)

            if (acc == "" || acc == "*") next

            # dedupe within cluster
            key = cluster "_" acc
            if (key in seen) next
            seen[key] = 1

            if (cluster in members) {
                members[cluster] = members[cluster] OFS acc
            } else {
                members[cluster] = acc
                order[++n] = cluster
            }
        }

        END {
            for (i = 1; i <= n; i++) {
                c = order[i]
                print "cluster_" c, members[c]
            }
        }
    ' "${uc}" \\
    > "${prefix}.cluster_accessions.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 || true)
    END_VERSIONS
    """
}
