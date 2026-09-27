nextflow.enable.dsl = 2

// Two real obipcr records: one with no primer mismatches, one with two. At a
// budget of 0 the second is discarded, so the smoke test shows the filter work.
params.obipcr_mismatches = 0

include { FILTER_OBIPCR } from '../../modules/filter_obipcr/main.nf'

workflow {
    def meta = [primer: 'ITS2', db: 'test_db']
    def raw_fasta = file("${projectDir}/dummy_raw.fasta")
    def primer_input = file("${projectDir}/primers")

    FILTER_OBIPCR(tuple(meta, raw_fasta, primer_input))

    FILTER_OBIPCR.out.stats.view()
}
