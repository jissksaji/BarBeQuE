nextflow.enable.dsl = 2

include { STAGE_FILE } from '../../modules/helper/stage_file/main.nf'

workflow {
    def f = file("${projectDir}/dummy_input.txt")
    
    STAGE_FILE(f)
    
    STAGE_FILE.out.staged_file.view()
}
