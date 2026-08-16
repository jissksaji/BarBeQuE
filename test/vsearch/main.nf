nextflow.enable.dsl = 2

include { VSEARCH_CLUSTER_FAST } from '../../modules/vsearch/cluster_fast/main.nf'

workflow {
    def meta = [primer: 'test_primer', db: 'test_db']
    def fa = file("${projectDir}/dummy.fasta")
    
    VSEARCH_CLUSTER_FAST(tuple(meta, fa))
    
    VSEARCH_CLUSTER_FAST.out.uc.view()
}
