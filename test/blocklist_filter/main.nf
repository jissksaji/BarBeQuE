nextflow.enable.dsl = 2

include { BLOCKLIST_FILTER } from '../../modules/helper/blocklist_filter/main.nf'

workflow {
    databases = channel.of(
        tuple([id: 'test_db'], file("${projectDir}/database.fasta"))
    )
    accession_taxonomy = channel.value(
        file("${projectDir}/accession_taxid.tsv")
    )
    blocklist = channel.value(file("${projectDir}/blocklist.txt"))

    BLOCKLIST_FILTER(databases, accession_taxonomy, blocklist)

    BLOCKLIST_FILTER.out.fasta.view()
    BLOCKLIST_FILTER.out.summary.view()
}
