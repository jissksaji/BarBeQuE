nextflow.enable.dsl = 2

params.outdir = "${projectDir}/results"

include { STREAMLIT } from '../../modules/streamlit/main.nf'

workflow {
    def results_dir = "${projectDir}/results"
    
    // create results dir
    file(params.outdir).mkdir()
    
    STREAMLIT(results_dir)
}
