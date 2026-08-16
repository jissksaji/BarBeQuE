nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../modules/blast/makeblastdb/main.nf'
include { BLASTN_CORE_NT } from '../../modules/blast/blastn_core_nt/main.nf'

workflow {
    def meta = [id: 'dummy_db']
    def fasta = file("${projectDir}/dummy_db.fasta")
    def taxid = file("${projectDir}/dummy_taxid.txt")
    
    BLAST_MAKEBLASTDB(tuple(meta, fasta, taxid))
    
    // Now test blastn_core_nt
    def query_meta = [id: 'query_1']
    // use the same fasta as query
    def query_fasta = file("${projectDir}/dummy_db.fasta")
    def taxid_val = 9606
    
    // wait for db to be created, pass it as val to blastn
    // BLAST_MAKEBLASTDB publishes the db dir; the actual blastdb prefix inside it
    // is the fasta filename (makeblastdb has no -out, so it defaults to -in's name)
    BLASTN_CORE_NT(tuple(query_meta, query_fasta), BLAST_MAKEBLASTDB.out.db.map{ db -> "${db[1]}/${fasta.name}" }, taxid_val)
    
    BLAST_MAKEBLASTDB.out.db.view()
    BLASTN_CORE_NT.out.tsv.view()
}
