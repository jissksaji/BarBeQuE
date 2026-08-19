nextflow.enable.dsl = 2

include { PARSE_PRIMERS } from '../../modules/parse_primers/main.nf'

workflow {
    // primers/ holds a plain single-pair FASTA, a file whose same-length variants collapse,
    // and a file holding two prefixes that must stay separate.
    PARSE_PRIMERS(
        channel.value([
            [id: 'test_primers', min: 100, max: 500],
            file("${projectDir}/primers"),
        ])
    )

    PARSE_PRIMERS.out.samplesheet.map { _meta, tsv -> tsv.text }.view()
}
