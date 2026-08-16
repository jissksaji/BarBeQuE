nextflow.enable.dsl = 2

params.obipcr_mismatches = 2

include { OBIPCR_INSILICOPCR } from '../../modules/obipcr/main.nf'

workflow {
    def meta = [
        primer: 'test_primer',
        db: 'dummy_db',
        fwd: 'GATC',
        rev: 'GATC',
        min: 10,
        max: 50
    ]
    def db_file = file("${projectDir}/dummy_db.fasta")
    
    OBIPCR_INSILICOPCR(tuple(meta, db_file))
    
    OBIPCR_INSILICOPCR.out.fasta.view()
}
