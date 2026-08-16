nextflow.enable.dsl = 2

include { DOWNLOAD_AND_COMBINE_DB } from '../../modules/download/main.nf'

workflow {
    def meta = [
        db: 'test_db',
        filename: 'combined_db.fasta',
        format: 'fasta',
        release: 'test',
    ]
    // Local fixture keeps this module test deterministic and usable offline.
    def urls = [file("${projectDir}/test.fasta").toUri().toString()]
    
    DOWNLOAD_AND_COMBINE_DB(channel.of(tuple(meta, urls)))
    
    DOWNLOAD_AND_COMBINE_DB.out.db.view()
}
