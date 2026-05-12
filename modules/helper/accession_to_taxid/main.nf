process CLUSTER_ACCESSIONS_TO_TAXID {

    tag "${meta.primer}|${meta.db}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cluster_accessions)
    path accession2taxid

    output:
    tuple val(meta), path("*.cluster_accession_taxid.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    set -euo pipefail

    read_taxmap() {
        if [[ "${accession2taxid}" == *.gz ]]; then
            gzip -cd "${accession2taxid}"
        else
            cat "${accession2taxid}"
        fi
    }

    long_clusters() {
        gawk -F'\\t' '
        BEGIN { OFS="\\t" }
        {
            cluster = \$1

            for (i = 2; i <= NF; i++) {
                acc = \$i
                gsub(/^ +| +\$/, "", acc)

                if (acc == "" || acc == "*") next

                lookup = acc
                sub(/\\.[0-9]+\$/, "", lookup)

                print cluster, acc, lookup
            }
        }
        ' "${cluster_accessions}"
    }

    gawk -F'\\t' '
    BEGIN {
        OFS="\\t"
        print "cluster_id", "accession", "taxid"
    }

    # File 1: wanted accessions
    ARGIND == 1 {
        want[\$1] = 1
        next
    }

    # File 2: accession2taxid
    ARGIND == 2 {
        if (FNR == 1 && \$1 == "accession") next
        if (NF < 2) next

        # NCBI format:
        # accession    accession.version    taxid    gi
        if (NF >= 3 && \$3 ~ /^[0-9]+\$/) {
            acc_nover = \$1
            acc_ver   = \$2
            tax        = \$3

            if (acc_nover in want) taxid[acc_nover] = tax
            if (acc_ver   in want) taxid[acc_ver]   = tax

            sub(/\\.[0-9]+\$/, "", acc_ver)
            if (acc_ver in want) taxid[acc_ver] = tax

            next
        }

        # Simple format:
        # accession    taxid
        if (NF >= 2 && \$2 ~ /^[0-9]+\$/) {
            acc = \$1
            tax = \$2

            if (acc in want) taxid[acc] = tax

            sub(/\\.[0-9]+\$/, "", acc)
            if (acc in want) taxid[acc] = tax

            next
        }
    }

    # File 3: long cluster table
    ARGIND == 3 {
        cluster = \$1
        acc     = \$2
        lookup  = \$3

        print cluster, acc, ((lookup in taxid) ? taxid[lookup] : "NA")
        next
    }
    ' \\
    <(long_clusters | cut -f3 | sort -u) \\
    <(read_taxmap) \\
    <(long_clusters) \\
    > "${prefix}.cluster_accession_taxid.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version 2>&1 | head -1 || true)
        gzip: \$(gzip --version 2>&1 | head -1 || true)
        sort: \$(sort --version 2>&1 | head -1 || true)
        cut: \$(cut --version 2>&1 | head -1 || true)
    END_VERSIONS
    """
}
