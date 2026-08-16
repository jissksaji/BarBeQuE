//
// This file holds functions to validate user-supplied arguments
//

import groovy.json.JsonOutput

class WorkflowPipeline {

    //
    // Check and validate parameters
    //
    public static void initialise( params, log) {
        if (params.build_references && !params.reference_base) {
            log.error '--build_references requires --reference_base'
            System.exit(1)
        }
         if (params.list_dbs) {
            println('Available databases:')
            println('===========================')
            params.references.databases.keySet().each { db ->
                def info = params.references.databases[db].description
                def cached = params.reference_base &&
                    new File("${params.reference_base}/barbeque/${params.reference_version}/filtered/${db}/${db}.cleaned.fasta").exists()
                println("Name: ${db}\tSource: ${info}\tFiltered cache: ${cached ? 'yes' : 'no'}")
                println('---------------------------')
            }
            System.exit(1)
        }
        if (params.list_primers) {
            def catalog
            try {
                catalog = PrimerCatalog.fetchCatalog()
            }
            catch (Exception e) {
                log.error("Could not fetch the primer-set catalog from bio-raum/FooDMe2: ${e.message}")
                System.exit(1)
            }
            println(JsonOutput.prettyPrint(JsonOutput.toJson(catalog)))
            System.exit(1)
        }
        if (params.hierarchical_clustering) {
            if (!params.custom_db) {
                log.info 'Standalone hierarchical clustering requires --custom_db <fasta>'
                System.exit(1)
            }
            if (params.dbs || params.input || params.primer_set) {
                log.info 'Standalone hierarchical clustering uses --custom_db directly; do not provide --dbs, --input, or --primer_set'
                System.exit(1)
            }
            if (!params.reference_base && (!params.taxdump || !params.accession_taxonomy)) {
                log.info 'Standalone hierarchical clustering needs taxonomy inputs: provide --reference_base, or provide both --taxdump and --accession_taxonomy'
                System.exit(1)
            }
            return
        }
        if (params.dbs && !params.run_name) {
            log.info 'Must provide a run_name (--run_name)'
            System.exit(1)
        }
        if (params.input && params.primer_set) {
            log.info 'Provide either --input or --primer_set, not both'
            System.exit(1)
        }
        if (!params.input && !params.primer_set && !params.build_references) {
            log.info "Pipeline requires a sample sheet / primer FASTA directory (--input) or a named primer set (--primer_set) as input"
            System.exit(1)
        }
        if (params.input) {
            def input_path = new File(params.input.toString())
            if (input_path.isDirectory() && (!params.primer_min || !params.primer_max)) {
                log.info "When --input is a primer FASTA directory, provide global amplicon bounds with --primer_min and --primer_max"
                System.exit(1)
            }
        }

    }

}
