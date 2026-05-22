process BUILD_DB_TAXIDS {

    tag "${meta.db}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path genbank2taxid

    output:
    tuple val(meta), path("*.db_taxids.tsv"), emit: taxids
    tuple val(meta), path("*.missing_accessions.tsv"), emit: missing
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.db}"
    """
    set -euo pipefail

    awk -F'\\t' '

        # PASS 1: FASTA — extract accession from header, load into lookup table
        /^>/ {
            accession = substr(\$1, 2)             # remove leading ">"
            sub(/\\.[0-9]+\$/, "", accession)      # strip version e.g. MK123456.1 -> MK123456
            db_accessions[accession] = 1
            next
        }

        # PASS 2: genbank2taxid — if accession matched, emit taxid and free from table
        {
            accession = \$1                        # col 1 = accession
            sub(/\\.[0-9]+\$/, "", accession)      # strip version to match PASS 1
            if (accession in db_accessions) {
                print \$2                          # col 2 = taxid
                delete db_accessions[accession]   # free from table to save memory
            }
        }

        # anything still in table was never found in genbank2taxid
        END {
            for (unmatched in db_accessions) {
                print unmatched > "${prefix}.missing_accessions.tsv"
            }
        }

    ' ${fasta} ${genbank2taxid} \\
    | sort -n -u \\
    > "${prefix}.db_taxids.tsv"

    touch "${prefix}.missing_accessions.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
