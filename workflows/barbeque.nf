// Modules
include { INPUT_CHECK } from './../modules/input_check'
include { MULTIQC } from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'
include { CUSTOM_DB_FILTER } from './../modules/seqkit/custom_db_filter'
//include { CRABS_INSILICOPCR }           from './../modules/crabs/insilico_pcr'
//include { CRABS_DEREPLICATE }           from './../modules/crabs/dereplicate'
//include { CRABS_FILTER }                from './../modules/crabs/filter'
//include { CRABS_SUBSET }                from './../modules/crabs/subset'
//include { CRABS_DIVERSITY_FIGURE }      from './../modules/crabs/diversity_figure'
include { VSEARCH_CLUSTER_FAST } from './../modules/vsearch/cluster_fast'
//include { CRABS_AMPLIFICATION_EFFICENCY_FIGURE } from './../modules/crabs/amplification_efficency_figure'
//include { CRABS_AMPLICON_LENGTH_FIGURE }from './../modules/crabs/amplicon_length_figure'
//include { HELPER_CLUSTER_CONSENSUS }    from './../modules/helper/cluster_consensus'
include { STAGE_FILE as STAGE_SAMPLESHEET } from './../modules/helper/stage_file'
//include { HELPER_CONSENSUS_HISTOGRAM }  from './../modules/helper/consensus_histogram'
//include { HELPER_TAXONOMIC_COVERAGE }   from './../modules/helper/taxonomic_coverage'
//include { HELPER_CONSENSUS_DISTRIBUTION } from './../modules/helper/consensus_distribution'
include { COMPUTE_BUFFER } from './../modules/seqkit/compute_buffer'
include { CUTADAPT_INSILICOPCR } from './../modules/cutadapt'
include { VSEARCH_DEREPLICATION } from './../modules/vsearch/dereplication'


workflow BARBEQUE {

    main:

    ch_multiqc_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : channel.value([])
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : channel.value([])

    ch_versions = channel.from([])
    multiqc_files = channel.from([])

    samplesheet = params.input ? channel.fromPath(file(params.input, checkIfExists: true)) : channel.value([])

    // The pre-installed taxdump folder

    // ch_taxdump = file(params.references.taxdump)

    //    pipeline_settings = channel.fromPath(dumpParametersToJSON(params.outdir)).collect()
    //
    //    // Check if the specified taxon is valid
    //    if (params.taxon) {
    //        taxon_valid = valid_taxon(params.taxon)
    //        if (!taxon_valid) {
    //          log.warn "Specified what appears to be an invalid taxon name - aborting!"
    //          System.exit(1)
    //        }
    //    }

    // the database to use - either pre-installed or user-provided
    // Pre-installed can be a list, coma-separated:  db1,db2,db3
    ch_dbs = channel.from([])
    these_dbs = []
    if (params.custom_db) {
        these_dbs << [["id": "custom"], file(params.custom_db, checkIfExists: true)]
    }
    else if (params.dbs) {
        valid_databases = params.references.databases.keySet()
        params.dbs
            .split(",")
            .collect { it.toLowerCase() }
            .each { db ->
                if (!valid_databases.contains(db)) {
                    log.info("Not a valid database: ${db}\nValid options are: ${valid_databases}\n")
                    System.exit(1)
                }
                these_dbs << [
                    [
                        "id": db
                    ],
                    file(params.references.databases[db].db, checkIfExists: true),
                ]
            }
    }
    ch_dbs = channel.fromList(these_dbs)

    ch_dbs.view { ">>> [1] DB: ${it}" }

    // Check if the samplesheet is valid
    INPUT_CHECK(samplesheet)

    INPUT_CHECK.out.primers.view { ">>> [2] PRIMER: ${it}" }

    // Copy the samplesheet to the results folder
    STAGE_SAMPLESHEET(samplesheet)

    //filter custom DB if provided
    if (params.custom_db) {
        CUSTOM_DB_FILTER(ch_dbs)
        ch_versions = ch_versions.mix(CUSTOM_DB_FILTER.out.versions)
        ch_dbs = CUSTOM_DB_FILTER.out.fasta

        ch_dbs.view { ">>> [3] FILTERED DB: ${it}" }
    }

    COMPUTE_BUFFER(ch_dbs)

    COMPUTE_BUFFER.out.buffersize.view { ">>> [4] BUFFER: ${it}" }

    /*
     Combine each primer set with all requested databases
     [ meta, database_meta, database_path ]
    */
    INPUT_CHECK.out.primers
        .combine(ch_dbs)
        .map { m, n, d ->
            [
                [
                    primer: m.primer,
                    fwd: m.fwd,
                    rev: m.rev,
                    min: m.min,
                    max: m.max,
                    db: n.id,
                ],
                d,
            ]
        }
        .set { ch_primers_with_db }

    ch_primers_with_db.view { ">>> [5] PRIMER+DB: ${it}" }

    CUTADAPT_INSILICOPCR(
        ch_primers_with_db.combine(COMPUTE_BUFFER.out.buffersize).view { ">>> [6] INTO CUTADAPT: ${it}" }
    )
    ch_versions = ch_versions.mix(CUTADAPT_INSILICOPCR.out.versions)

    CUTADAPT_INSILICOPCR.out.fasta.view { ">>> [7] OUT CUTADAPT: ${it}" }

    CUTADAPT_INSILICOPCR.out.fasta
        .branch { m, f ->
            valid: file(f).size() > 0
            invalid: file(f).size() == 0
        }
        .set { ch_insilico_by_status }

    ch_insilico_by_status.valid.view { ">>> [8] VALID: ${it}" }

    ch_insilico_by_status.invalid.subscribe { m, t ->
        log.warn("${m.primer} did not produce any pcr products, stopping primer set")
    }
    VSEARCH_DEREPLICATION(ch_insilico_by_status.valid)
    ch_versions = ch_versions.mix(VSEARCH_DEREPLICATION.out.versions)

    VSEARCH_DEREPLICATION.out.fasta.view { ">>> [9] DEREPLICATION: ${it}" }

    VSEARCH_CLUSTER_FAST(VSEARCH_DEREPLICATION.out.fasta)

    //VSEARCH_CLUSTER_FAST.out.fasta.view { ">>> [10] FASTA CENTROIDS: ${it}" }
    VSEARCH_CLUSTER_FAST.out.uc.view { ">>> [10] UC CLUSTERING: ${it}" }


    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Combine by meta dict to generate separate reports for each primer-db combination
    multiqc_by_set = multiqc_files.groupTuple(by: 0)

    MULTIQC(
        multiqc_by_set,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
        ch_multiqc_config,
        ch_multiqc_logo,
    )

    emit:
    qc = MULTIQC.out.html
}

//
//    // perform insilico pcr, takes: [meta, database]
//    CRABS_INSILICOPCR(
//        ch_primers_with_db
//    )
//    ch_versions = ch_versions.mix(CRABS_INSILICOPCR.out.versions)
//
//    CRABS_INSILICOPCR.out.txt.branch { m,t ->
//        valid: file(t).size() > 0
//        invalid: file(t).size() == 0
//    }.set { ch_insilico_by_status }
//
//    ch_insilico_by_status.invalid.subscribe { m,t ->
//        log.warn "${m.primer} did not produce any pcr products, stopping primer set"
//    }
//
//    // dereplicate in-silico amplicons, takes [meta, txt]
//    CRABS_DEREPLICATE(
//        ch_insilico_by_status.valid
//    )
//    ch_versions = ch_versions.mix(CRABS_DEREPLICATE.out.versions)
//
//    // Filter hits, takes [meta, txt]
//    CRABS_FILTER(
//        CRABS_DEREPLICATE.out.txt
//    )
//    ch_versions = ch_versions.mix(CRABS_FILTER.out.versions)
//
//    // fast clustering of crabs OTUs
//    VSEARCH_CLUSTER_FAST(
//        CRABS_FILTER.out.fasta
//    )
//    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_FAST.out.versions)
//
//    // Cluster consensus
//    HELPER_CLUSTER_CONSENSUS(
//        VSEARCH_CLUSTER_FAST.out.uc.join(CRABS_FILTER.out.txt),
//        ch_taxdump
//    )
//    ch_versions = ch_versions.mix(HELPER_CLUSTER_CONSENSUS.out.versions)
//
//    HELPER_CLUSTER_CONSENSUS.out.txt.map { m, t ->
//        tuple(m.db, m, t)
//    }.combine(
//        ch_dbs.map { n, d ->
//            tuple(n.id, d)
//        }, by: 0
//    ).map { k, m, t, d ->
//        tuple(m, t, d)
//    }.set { ch_cluster_with_db }
//
//    // Amplicon size distribution
//    HELPER_CONSENSUS_DISTRIBUTION(
//        ch_cluster_with_db
//    )
//    ch_versions = ch_versions.mix(HELPER_CONSENSUS_DISTRIBUTION.out.versions)
//    multiqc_files = multiqc_files.mix(HELPER_CONSENSUS_DISTRIBUTION.out.json)
//
//    // convert the consensus file into a histogram of amplicon lengths
//    HELPER_CONSENSUS_HISTOGRAM(
//        HELPER_CLUSTER_CONSENSUS.out.txt
//    )
//    multiqc_files = multiqc_files.mix(HELPER_CONSENSUS_HISTOGRAM.out.json)
//    ch_versions = ch_versions.mix(HELPER_CONSENSUS_HISTOGRAM.out.versions)
//
//    // If a taxon is provided, perform additional visualisation/filtering
//    if (params.taxon) {
//
//        // Analyse the coverage of the desired taxonomic level
//        HELPER_TAXONOMIC_COVERAGE(
//            HELPER_CLUSTER_CONSENSUS.out.txt.map { m,t ->
//                tuple(m.db, m, t)
//            }.combine(
//                ch_dbs.map { m,d ->
//                    tuple(m.id,d)
//                }, by: 0
//            ).map { k, m, s, d ->
//                tuple(m,s,d)
//            },
//            params.taxon
//        )
//
//        // Generate a subset based on the --taxon argument
//        CRABS_SUBSET(
//            CRABS_FILTER.out.txt,
//            params.taxon
//        )
//
//        CRABS_SUBSET.out.txt.branch { m,t ->
//            valid: t.size() > 0
//            invalid: t.size() == 0
//        }.set { ch_subset_by_status }
//
//        ch_subset_by_status.invalid.subscribe {m,t ->
//            log.warn "No hits left after subsetting ${m.primer} with ${params.taxon} - stopping."
//        }
//        
//        // Visualize the length distribution of putative amplicons
//        CRABS_AMPLICON_LENGTH_FIGURE(
//            ch_subset_by_status.valid
//        )
//
//        // Visualize diversity of amplicons
//        CRABS_DIVERSITY_FIGURE(
//            ch_subset_by_status.valid
//        )
//
//        // Combine each subset with the correct database
//        ch_subset_by_status.valid.map {m, s ->
//            tuple(m.db,m,s)
//        }.combine(
//            ch_dbs.map { m,d ->
//                tuple(m.id,d)
//            }, by: 0
//        ).map { k, m, s, d ->
//            tuple(m,s,d)
//        }.set { ch_amplicons_with_db }
//
//        // visualize amplification efficency
//        CRABS_AMPLIFICATION_EFFICENCY_FIGURE(
//            ch_amplicons_with_db,
//            params.taxon
//        )
//    }
//
//    CUSTOM_DUMPSOFTWAREVERSIONS(
//        ch_versions.unique().collectFile(name: 'collated_versions.yml')
//    )
//
//    // Combine by meta dict to generate separate reports for each primer-db combination
//    multiqc_by_set = multiqc_files.groupTuple(by: 0)
//
//    MULTIQC(
//        multiqc_by_set,
//        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
//        ch_multiqc_config,
//        ch_multiqc_logo
//    )
//
//    emit:
//    qc = MULTIQC.out.html
//}
//

//// turn the params map to a JSON file
def dumpParametersToJSON(outdir) {
    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
    def filename = "params_${timestamp}.json"
    def temp_pf = new File(workflow.launchDir.toString(), ".${filename}")
    def jsonStr = groovy.json.JsonOutput.toJson(params)
    temp_pf.text = groovy.json.JsonOutput.prettyPrint(jsonStr)

    nextflow.extension.FilesEx.copyTo(temp_pf.toPath(), "${outdir}/pipeline_info/params_${timestamp}.json")
    temp_pf.delete()
    return file("${outdir}/pipeline_info/params_${timestamp}.json")
}
