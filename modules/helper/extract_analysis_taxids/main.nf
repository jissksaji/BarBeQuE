process EXTRACT_ANALYSIS_TAXIDS {

    publishDir "${params.outdir}/analysis_taxids", mode: 'copy'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta), path(accession_taxid)
    path taxdump
    val taxids

    output:
    path "*_*.fasta"

    script:
    """
    for taxid in ${taxids}; do

        # Get the scientific name of the TaxID
        # Example: 4442 becomes Camellia_sinensis
        # Replace spaces with underscores for the filename

        name=\$(echo \$taxid | taxonkit lineage -n -L --data-dir ${taxdump} | cut -f2 | tr ' ' '_')

        # Get the taxonomic rank of the TaxID
        # Example: species, genus, family

        rank=\$(echo \$taxid | taxonkit lineage -r -L --data-dir ${taxdump} | cut -f2)

        # Find the TaxID and all TaxIDs below it in the taxonomy tree
        # Example: a genus includes all species belonging to that genus
        # Save all these TaxIDs in taxids.txt

        taxonkit list --ids \$taxid --indent "" --data-dir ${taxdump} > taxids.txt

        # Find accessions associated with the TaxIDs in taxids.txt
        # grep searches for matching TaxIDs in the accession-taxid file
        # cut -f1 keeps only the first column (accession numbers)
        # Save the matching accessions in accessions.txt

        grep -wFf taxids.txt ${accession_taxid} | cut -f1 > accessions.txt

        # Extract FASTA sequences whose IDs match accessions.txt
        # --id-regexp extracts the accession without the version number
        # Example: AB123456.1 becomes AB123456 for matching
        # Keep the original FASTA headers and complete sequences
        # Save the sequences in a separate FASTA file for each TaxID
        
        seqkit grep --id-regexp '^(\\w+)' -f accessions.txt ${fasta} > \${taxid}_\${name}_\${rank}.fasta

    done
    """
}
