nextflow.enable.dsl = 2

include { FASTP } from '../../modules/fastp/main.nf'

workflow {
    def meta = [id: 'test_sample', single_end: false]
    def reads = tuple(
        file("${projectDir}/test_R1.fastq"),
        file("${projectDir}/test_R2.fastq")
    )
    
    FASTP(tuple(meta, reads))
    
    FASTP.out.reads.view()
}
