nextflow.enable.dsl = 2

include { PARSE_OBIPCR } from '../../modules/parse_obipcr/main.nf'

workflow {
    def meta = [primer: 'test_primer', db: 'test_db']
    def fasta = file("${projectDir}/dummy.fasta")
    
    PARSE_OBIPCR(tuple(meta, fasta))
    
    PARSE_OBIPCR.out.tsv.view()
}
