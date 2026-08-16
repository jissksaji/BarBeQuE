nextflow.enable.dsl = 2

params.outdir = "${projectDir}/results"

include { SAMTOOLS_FAIDX } from '../../modules/samtools/faidx/main.nf'

workflow {
    def meta = [id: 'test_samtools']
    def fasta = file("${projectDir}/dummy.fasta")
    
    SAMTOOLS_FAIDX(tuple(meta, fasta))
    
    SAMTOOLS_FAIDX.out.fai.view()
}
