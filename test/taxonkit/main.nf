nextflow.enable.dsl = 2

include { TAXONKIT_LCA } from '../../modules/taxonkit/lca/main.nf'

workflow {
    def meta = [primer: 'test_primer', db: 'test_db']
    def cluster_taxids = file("${projectDir}/dummy_cluster_taxids.tsv")
    def taxdump = file("${projectDir}/../db_distribution/taxdump")
    
    TAXONKIT_LCA(tuple(meta, cluster_taxids), taxdump)
    
    TAXONKIT_LCA.out.lca.view()
}
