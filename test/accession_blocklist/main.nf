nextflow.enable.dsl = 2

include { ACCESSION_BLOCKLIST } from '../../modules/helper/accession_blocklist/main.nf'

workflow {
    def meta = [primer: 'test_primer', db: 'test_db']
    amplicons = channel.value(tuple(
        meta,
        file("${projectDir}/amplicons.fasta"),
        file("${projectDir}/parsed.tsv"),
    ))
    accession_blocklist = channel.value(
        file("${projectDir}/accession_blocklist.txt")
    )

    ACCESSION_BLOCKLIST(amplicons, accession_blocklist)

    ACCESSION_BLOCKLIST.out.fasta.view()
    ACCESSION_BLOCKLIST.out.tsv.view()
    ACCESSION_BLOCKLIST.out.summary.view()
}
