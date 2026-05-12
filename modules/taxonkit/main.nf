process ADD_TAXONOMY_TAXONKIT {

    tag "${meta.primer}|${meta.db}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cluster_accession_taxid)
    path taxdump_dir

    output:
    tuple val(meta), path("*.cluster_accession_taxonomy.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    """
    set -euo pipefail

    # 1. Extract unique taxids from cluster/accession/taxid table
    awk -F '\\t' '
        NR > 1 && \$3 != "NA" && \$3 ~ /^[0-9]+\$/ {
            print \$3
        }
    ' "${cluster_accession_taxid}" | sort -u > taxids.txt

    # 2. Convert taxids to full + ranked lineage
    taxonkit lineage --data-dir "${taxdump_dir}" taxids.txt \\
    | taxonkit reformat \\
        --data-dir "${taxdump_dir}" \\
        -r NA \\
        -f "{k}\\t{p}\\t{c}\\t{o}\\t{f}\\t{g}\\t{s}" \\
    > taxid_taxonomy.tsv

    # taxid_taxonomy.tsv columns:
    # taxid full_lineage superkingdom phylum class order family genus species

    # 3. Join taxonomy back to original table
    awk -F '\\t' '
    BEGIN { OFS="\\t" }

    NR == FNR {
        tax[\$1] = \$2 OFS \$3 OFS \$4 OFS \$5 OFS \$6 OFS \$7 OFS \$8 OFS \$9
        next
    }

    FNR == 1 {
        print \$0, "full_lineage", "superkingdom", "phylum", "class", "order", "family", "genus", "species"
        next
    }

    {
        taxid = \$3
        print \$0, ((taxid in tax) ? tax[taxid] : "NA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA")
    }
    ' taxid_taxonomy.tsv "${cluster_accession_taxid}" \\
    > "${prefix}.cluster_accession_taxonomy.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$(taxonkit version 2>&1 | head -1 || true)
        awk: \$(awk --version 2>&1 | head -1 || true)
        sort: \$(sort --version 2>&1 | head -1 || true)
    END_VERSIONS
    """
}
