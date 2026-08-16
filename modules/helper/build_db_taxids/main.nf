process BUILD_DB_TAXIDS {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"


    input:
    tuple val(meta), path(fasta)
    path genbank2taxid

    output:
    tuple val(meta), path("*.db_taxids.tsv"), emit: taxids
    tuple val(meta), path("*.db_taxids_counts.tsv"), emit: taxids_counts
    tuple val(meta), path("*.accession_taxid.tsv"), emit: accession_taxid
    tuple val(meta), path("*.missing_accessions.tsv"), emit: missing
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    # Pre-extract accessions to drastically speed up processing and avoid awk memory limits on huge fastas
    # cut -d';' keeps SINTAX-formatted headers (>ACC;tax=k:...,s:...;) working alongside plain
    # GenBank headers (>ACC.1 description)
    grep "^>" ${fasta} | cut -d' ' -f1 | cut -d'>' -f2 | cut -d';' -f1 > db_accessions.txt

    awk -F'\\t' -v acc_taxid_out="${prefix}.accession_taxid.tsv" '
        # PASS 1: DB accessions — load into lookup table
        NR == FNR {
            sub(/(\\.[0-9]+)+\$/, "", \$1)   # strip version/range e.g. MK123456.1 or AY846379.1.1791 -> base accession
            db_accessions[\$1] = 1
            next
        }

        # PASS 2: genbank2taxid — if accession matched, emit taxid and free from table
        {
            accession = \$1                        # col 1 = accession
            sub(/(\\.[0-9]+)+\$/, "", accession)   # strip version/range to match PASS 1
            if (accession in db_accessions) {
                taxid = (NF >= 3) ? \$3 : \$2
                print accession"\\t"taxid > acc_taxid_out   # so later per-primer joins can reuse this instead of re-scanning genbank2taxid
                print taxid                          # col 2 or 3 = taxid
                delete db_accessions[accession]   # free from table to save memory
            }
        }

        # anything still in table was never found in genbank2taxid
        END {
            for (unmatched in db_accessions) {
                print unmatched > "${prefix}.missing_accessions.tsv"
            }
        }
    ' db_accessions.txt ${genbank2taxid} \\
    | sort -n \\
    | tee >(sort -u > "${prefix}.db_taxids.tsv") \\
    | uniq -c \\
    | awk '{print \$2"\\t"\$1}' \\
    > "${prefix}.db_taxids_counts.tsv"

    rm db_accessions.txt

    touch "${prefix}.missing_accessions.tsv"
    touch "${prefix}.accession_taxid.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n 1 || echo "unknown")
    END_VERSIONS
    """
}
