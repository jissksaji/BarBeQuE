process JOIN_ACCESSION_TAXONOMY {

    tag "${meta.primer}|${meta.db}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cluster_accessions)
    path master_lookup

    output:
    tuple val(meta), path("*.cluster_taxonomy.tsv"), emit: tsv
    tuple val(meta), path("*.missing_accessions.tsv"), emit: missing
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"
    """
set -euo pipefail

awk -F'\\t' '
    BEGIN { OFS = "\\t" }

    # First file: load cluster_id + accession into hash
    NR == FNR {
        cluster[\$2] = \$1
        wanted[\$2] = 1
        next
    }

    # Second file: data rows — join if accession is wanted
    (\$1 in cluster) {
        print cluster[\$1], \$0
        found[\$1] = 1
    }

    END {
        # Write missing accessions to separate file
        for (acc in wanted) {
            if (!(acc in found)) {
                print cluster[acc] "\\t" acc > "${prefix}.missing_accessions.tsv"
            }
        }
    }
' ${cluster_accessions} ${master_lookup} \\
> "${prefix}.cluster_taxonomy.tsv"

# Ensure missing file exists even if empty
touch "${prefix}.missing_accessions.tsv"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    awk: \$(awk --version 2>&1 | head -n 1 || echo "unknown")
END_VERSIONS
"""
}
