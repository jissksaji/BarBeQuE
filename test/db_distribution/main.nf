nextflow.enable.dsl = 2

include { DB_DISTRIBUTION } from '../../modules/db_distribution/main.nf'

workflow {
    def meta = [id: 'test_db_distribution']
    def taxids_counts = file("${projectDir}/dummy_counts.tsv")
    def taxdump = file("${projectDir}/taxdump")
    
    DB_DISTRIBUTION(tuple(meta, taxids_counts), taxdump)
    
    DB_DISTRIBUTION.out.distribution.view()
}
