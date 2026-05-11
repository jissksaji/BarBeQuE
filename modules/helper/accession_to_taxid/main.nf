process CLUSTER_ACCESSIONS_TO_TAXID {

    tag "${meta.primer}|${meta.db}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cluster_accessions)
    path accession2taxid

    output:
    tuple val(meta), path("*.accession_taxid.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    set -euo pipefail

    # Extract unique accessions from cluster file
    awk -F'\\t' '
        {
            for (i = 2; i <= NF; i++) {
                acc = \$i
                sub(/\\.[0-9]+\$/, "", acc)
                if (acc != "" && acc != "*") print acc
            }
        }
    ' "${cluster_accessions}" \\
    | sort -u \\
    > wanted.txt

    # Lookup taxids
    if [[ "${accession2taxid}" == *.gz ]]; then
        READ_CMD="gzip -cd"
    else
        READ_CMD="cat"
    fi

    \$READ_CMD "${accession2taxid}" \\
    | awk -F'\\t' '
        NR == FNR { want[\$1] = 1; next }
        FNR == 1 && \$1 == "accession" { next }
        NF < 2 { next }
        {
            acc = \$1
            sub(/\\.[0-9]+\$/, "", acc)
            if (!(acc in want)) next

            if (NF >= 3 && \$3 ~ /^[0-9]+\$/) print acc "\\t" \$3
            else if (\$2 ~ /^[0-9]+\$/) print acc "\\t" \$2
        }
    ' wanted.txt - \\
    | sort -u \\
    > "${prefix}.accession_taxid.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 || true)
        gzip: \$(gzip --version 2>&1 | head -1 || true)
        sort: \$(sort --version 2>&1 | head -1 || true)
    END_VERSIONS
    """
}
