nextflow.enable.dsl = 2

include { PRIMER_DISAMBIGUATE } from '../../modules/primer_disambiguate/main.nf'

workflow {
    def fa = file("${projectDir}/dummy_primer.fasta")
    
    PRIMER_DISAMBIGUATE(fa)
    
    PRIMER_DISAMBIGUATE.out.fasta.view()
}
