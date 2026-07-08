/*
Import modules
*/
include { INPUT_CHECK } from './../modules/input_check'
include { MULTIQC } from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'
include { STAGE_FILE as STAGE_SAMPLESHEET } from './../modules/helper/stage_file'

/*
Helper Modules
*/
include { COMPUTE_BUFFER } from './../modules/seqkit/compute_buffer'
include { BUILD_DB_TAXIDS } from './../modules/helper/build_db_taxids'
include { DB_DISTRIBUTION } from './../modules/db_distribution/main.nf'
include { AMPLICON_LENGTH } from './../modules/seqkit/amplicon_lengths'
include { TAXONOMIC_COVERAGE_PLOT } from './../modules/helper/taxonomic_coverage_plot/main'
include { COMPLETENESS_TABLE } from './../modules/helper/completeness_table/main'
include { SPECIES_REPRESENTATION } from './../modules/helper/species_representation/main'
include { MASK } from './../modules/mask/main.nf'
include { PARSE_UC } from './../modules/helper/parse_uc'
include { JOIN_ACCESSION_TAXONOMY } from './../modules/helper/join_accession_taxonomy'
/*
Core Modules
*/
include { CUTADAPT_INSILICOPCR } from './../modules/cutadapt'
include { OBIPCR_INSILICOPCR } from './../modules/obipcr'
include { PARSE_OBIPCR } from './../modules/parse_obipcr/main.nf'
include { VSEARCH_DEREPLICATION } from './../modules/vsearch/dereplication'
include { VSEARCH_CLUSTER_FAST } from './../modules/vsearch/cluster_fast'
include { TAXONOMIC_COVERAGE } from './../modules/helper/taxonomic_coverage'
include { CLUSTER_CONSENSUS } from './../modules/helper/cluster_consensus'



workflow BARBEQUE {

    take:
    ch_dbs
    ch_db_versions

    main:
    samplesheet = params.input ? channel.fromPath(file(params.input, checkIfExists: true)) : channel.value([])
    ch_multiqc_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : channel.value([])
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : channel.value([])
    ch_versions = ch_db_versions
    multiqc_files = channel.from([])
    def taxdump_path = params.taxdump ?: params.references.taxdump
    ch_taxdump = channel.value(
        file(taxdump_path, checkIfExists: !params.build_references)
    )

    //channel for acession_to_taxonomy(basically genbank2taxid)
    ch_accession_taxonomy = channel.from([])

    def accession_file = params.accession_taxonomy
    if (!accession_file && taxdump_path) {
        accession_file = files("${taxdump_path}/*{accession2taxid,genbank2taxid,nucl_}*").find()
    }

    if (accession_file) {
        ch_accession_taxonomy = channel.value(
            file(accession_file, checkIfExists: true)
        )
    }
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

    //ch_dbs.view { ">>> [1] DB: ${it}" }


    // Check if the samplesheet is valid
    INPUT_CHECK(samplesheet)

    // Copy the samplesheet to the results folder
    STAGE_SAMPLESHEET(samplesheet)






    //COMPUTE_BUFFER.out.buffersize.view { ">>> [4] BUFFER: ${it}" }

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

    //ch_primers_with_db.view { ">>> [5] PRIMER+DB: ${it}" }

    def insilico_tool = (params.insilico_tool ?: 'obipcr').toLowerCase()
    if (!(insilico_tool in ['obipcr', 'cutadapt'])) {
        log.error("Invalid --insilico_tool '${insilico_tool}' - must be 'obipcr' or 'cutadapt'")
        System.exit(1)
    }

    if (insilico_tool == 'cutadapt') {
        COMPUTE_BUFFER(ch_dbs)
        CUTADAPT_INSILICOPCR(
            ch_primers_with_db.combine(COMPUTE_BUFFER.out.buffersize)
        )
        ch_versions = ch_versions.mix(CUTADAPT_INSILICOPCR.out.versions)
        ch_insilico_fasta = CUTADAPT_INSILICOPCR.out.fasta
    }
    else {
        OBIPCR_INSILICOPCR(
            ch_primers_with_db
        )
        ch_versions = ch_versions.mix(OBIPCR_INSILICOPCR.out.versions)
        ch_insilico_fasta = OBIPCR_INSILICOPCR.out.fasta
        
        PARSE_OBIPCR(
            OBIPCR_INSILICOPCR.out.raw_fasta
        )
        ch_versions = ch_versions.mix(PARSE_OBIPCR.out.versions)
    }

    ch_insilico_fasta
        .branch { _m, f ->
            valid: file(f).size() > 0
            invalid: file(f).size() == 0
        }
        .set { ch_insilico_by_status }

    if (params.mask) {
        MASK(ch_insilico_by_status.valid)
        ch_amplicons = MASK.out.fasta
    }
    else {
        ch_amplicons = ch_insilico_by_status.valid
    }

    AMPLICON_LENGTH(ch_amplicons)
    ch_versions = ch_versions.mix(AMPLICON_LENGTH.out.versions)

    VSEARCH_CLUSTER_FAST(ch_amplicons)

    //VSEARCH_CLUSTER_FAST.out.fasta.view { ">>> [10] FASTA CENTROIDS: ${it}" }
    //VSEARCH_CLUSTER_FAST.out.uc.view { ">>> [10] UC CLUSTERING: ${it}" }

    PARSE_UC(VSEARCH_CLUSTER_FAST.out.uc)
    ch_versions = ch_versions.mix(PARSE_UC.out.versions)

    // Extract the accession->taxid mapping once per db (not per primer) - this is the
    // expensive full scan of the (potentially huge) master accession2taxid/genbank2taxid file.
    BUILD_DB_TAXIDS(ch_dbs, ch_accession_taxonomy)
    ch_versions = ch_versions.mix(BUILD_DB_TAXIDS.out.versions)

    JOIN_ACCESSION_TAXONOMY(
        PARSE_UC.out.tsv
            .map { m, tsv -> tuple(m.db, m, tsv) }
            .combine(
                BUILD_DB_TAXIDS.out.accession_taxid.map { m, f -> tuple(m.id, f) },
                by: 0
            )
            .map { _db_id, m, tsv, f -> tuple(m, tsv, f) }
    )
    ch_versions = ch_versions.mix(JOIN_ACCESSION_TAXONOMY.out.versions)

    CLUSTER_CONSENSUS(
        JOIN_ACCESSION_TAXONOMY.out.tsv,
        ch_taxdump,
    )
    ch_versions = ch_versions.mix(CLUSTER_CONSENSUS.out.versions)

    DB_DISTRIBUTION(
        BUILD_DB_TAXIDS.out.taxids_counts,
        ch_taxdump
    )
    ch_versions = ch_versions.mix(DB_DISTRIBUTION.out.versions)

    if (params.completeness_table) {
        COMPLETENESS_TABLE(
            BUILD_DB_TAXIDS.out.taxids_counts
        )
        ch_versions = ch_versions.mix(COMPLETENESS_TABLE.out.versions)
    }

    if (params.taxon) {
        TAXONOMIC_COVERAGE(
            CLUSTER_CONSENSUS.out.tsv.map { m, tsv ->
                tuple(m.db, m, tsv)
            }.combine(
                BUILD_DB_TAXIDS.out.taxids.map { m, taxids ->
                    tuple(m.id, taxids)
                },
                by: 0
            ).map { db_id, meta, tsv, taxids ->
                tuple(meta, tsv, taxids)
            },
            params.taxon,
        )

        ch_versions = ch_versions.mix(TAXONOMIC_COVERAGE.out.versions)

        SPECIES_REPRESENTATION(
            TAXONOMIC_COVERAGE.out.tsv.join(CLUSTER_CONSENSUS.out.tsv, by: 0)
        )
        ch_versions = ch_versions.mix(SPECIES_REPRESENTATION.out.versions)
    }

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
    consensus = CLUSTER_CONSENSUS.out.tsv
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
//        // Visualize the length distribution of putati