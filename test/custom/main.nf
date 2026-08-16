nextflow.enable.dsl = 2

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/custom/dumpsoftwareversions/main.nf'

workflow {
    def versions_file = file("${projectDir}/dummy_versions.yml")
    
    CUSTOM_DUMPSOFTWAREVERSIONS(channel.of(versions_file).collect())
    
    CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.view()
}
