nextflow.enable.dsl = 2

params.read_length = 150
params.single_end = false

include { MASK } from '../../modules/mask/main.nf'

workflow {
    def meta = [primer: 'test_primer', db: 'test_db']
    def fasta = file("${projectDir}/dummy.fasta")
    
    MASK(tuple(meta, fasta))
    
    MASK.out.fasta.view()
}
