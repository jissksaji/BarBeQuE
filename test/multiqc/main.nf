nextflow.enable.dsl = 2

include { MULTIQC } from '../../modules/multiqc/main.nf'

workflow {
    def meta = [primer: 'test_primer', db: 'test_db']
    
    // We pass a json file simulating fastp output, config, logo, software
    def input_files = tuple(meta, [file("${projectDir}/dummy_fastp.json")])
    def software = file("${projectDir}/dummy_software.yml")
    def config = file("${projectDir}/dummy_config.yml")
    def logo = file("${projectDir}/dummy_logo.png")
    
    MULTIQC(input_files, software, config, logo)
    
    MULTIQC.out.html.view()
}
