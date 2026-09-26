include { PRIMER_SET } from './../subworkflows/primer_set/main'
include { PARSE_PRIMERS } from './../modules/parse_primers/main'
include { INPUT_CHECK } from './../modules/input_check'
include { STAGE_FILE as STAGE_SAMPLESHEET } from './../modules/helper/stage_file/main'
include { OBIPCR_INSILICOPCR } from './../modules/obipcr/main'
include { FILTER_OBIPCR } from './../modules/filter_obipcr/main'
include { PARSE_OBIPCR } from './../modules/parse_obipcr/main'
include { OBIPCR_REPORT } from './../modules/obipcr_report/main'
include { VSEARCH_DEREPLICATION } from './../modules/vsearch/dereplication/main'
include { MASK } from './../modules/mask/main'
include { BUILD_DB_TAXIDS } from './../modules/helper/build_db_taxids/main'
include { VSEARCH_CLUSTER_FAST } from './../modules/vsearch/cluster_fast/main'
include { PARSE_UC } from './../modules/helper/parse_uc/main'
include { JOIN_ACCESSION_TAXONOMY } from './../modules/helper/join_accession_taxonomy/main'
include { CLUSTER_CONSENSUS } from './../modules/helper/cluster_consensus/main'
include { CLUSTER_CONSENSUS_REPORTING } from './../modules/cluster_consensus_reporting/main'
include { DB_DISTRIBUTION } from './../modules/db_distribution/main'
include { TAXONOMIC_COVERAGE } from './../modules/helper/taxonomic_coverage/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions/main'
include { REPORTING } from './../subworkflows/reporting/main'

// Pair primer results with database lookups by database name.
def combine_by_db(ch_per_primer, ch_per_db) {
    return ch_per_primer
        .map { meta, f -> tuple(meta.db, meta, f) }
        .combine(ch_per_db.map { meta, lookup -> tuple(meta.id, lookup) }, by: 0)
        .map { _db_id, meta, f, lookup -> tuple(meta, f, lookup) }
}

workflow BARBEQUE {
    take:
    ch_dbs
    ch_db_versions
    ch_taxdump
    ch_accession_taxonomy

    main:
    ch_versions = ch_db_versions
    multiqc_files = channel.empty()

    // Set only for primer FASTA input, where the uncollapsed primers still exist
    // and FILTER_OBIPCR can check amplicons against them.
    primer_input = null

    // Load primers from the catalog, FASTA input, or a samplesheet.
    if (params.primer_set) {
        PRIMER_SET()
        samplesheet = PRIMER_SET.out.samplesheet
        ch_versions = ch_versions.mix(PRIMER_SET.out.versions)
    }
    else if (WorkflowPipeline.isFastaInput(params.input)) {
        primer_input = file(params.input, checkIfExists: true)

        // Convert primer FASTA input to a samplesheet.
        PARSE_PRIMERS(
            channel.value(
                [
                    [id: 'primers', min: params.primer_min, max: params.primer_max],
                    primer_input,
                ]
            )
        )
        ch_versions = ch_versions.mix(PARSE_PRIMERS.out.versions)
        samplesheet = PARSE_PRIMERS.out.samplesheet.map { _meta, tsv -> tsv }
    }
    else {
        samplesheet = channel.fromPath(file(params.input, checkIfExists: true))
    }

    // Validate the samplesheet and make primer names path-safe.
    INPUT_CHECK(samplesheet)

    // Copy the samplesheet to pipeline_info/ through publishDir.
    STAGE_SAMPLESHEET(samplesheet)

    ch_primers = INPUT_CHECK.out.primers

    // Create one in-silico PCR input for each primer and database pair.
    ch_primers
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


    OBIPCR_INSILICOPCR(
        ch_primers_with_db
    )
    ch_versions = ch_versions.mix(OBIPCR_INSILICOPCR.out.versions)

    // PARSE_PRIMERS combines primer variants into one degenerate primer, which also
    // matches base combinations no real primer has. With primer FASTA input the
    // originals are still around, so drop the amplicons only the combined primer
    // could explain. Both obipcr outputs are chosen together here, so the report
    // and the amplicons sent downstream always describe the same sequences.
    if (primer_input != null && WorkflowPipeline.enabled(params.filter_collapsed_primers)) {
        FILTER_OBIPCR(
            OBIPCR_INSILICOPCR.out.raw_fasta.map { meta, raw -> [meta, raw, primer_input] }
        )
        ch_versions = ch_versions.mix(FILTER_OBIPCR.out.versions)
        ch_obipcr_raw = FILTER_OBIPCR.out.raw_fasta
        ch_obipcr_fasta = FILTER_OBIPCR.out.fasta
    }
    else {
        ch_obipcr_raw = OBIPCR_INSILICOPCR.out.raw_fasta
        ch_obipcr_fasta = OBIPCR_INSILICOPCR.out.fasta
    }

    // Parse raw obipcr headers into a per-amplicon TSV.
    PARSE_OBIPCR(
        ch_obipcr_raw
    )
    ch_versions = ch_versions.mix(PARSE_OBIPCR.out.versions)

    OBIPCR_REPORT(PARSE_OBIPCR.out.tsv)
    ch_versions = ch_versions.mix(OBIPCR_REPORT.out.versions)
    // Emit each MultiQC JSON file separately.
    multiqc_files = multiqc_files.mix(OBIPCR_REPORT.out.mqc.transpose())

    // Keep the latest amplicon output as optional stages run.
    ch_amplicons = ch_obipcr_fasta

    // Warn about and remove empty amplicon results.
    ch_amplicons = ch_amplicons.filter { meta, fasta ->
        if (fasta.size() == 0) {
            log.warn("Primer produced no amplicons: ${meta.primer} (${meta.db})")
        }
        fasta.size() > 0
    }

    if (WorkflowPipeline.enabled(params.dereplicate_amplicons)) {
        VSEARCH_DEREPLICATION(ch_amplicons)
        ch_versions = ch_versions.mix(VSEARCH_DEREPLICATION.out.versions)
        ch_amplicons = VSEARCH_DEREPLICATION.out.fasta
    }

    if (WorkflowPipeline.enabled(params.mask)) {
        // Mask amplicons to simulate the configured read length.
        MASK(ch_amplicons)
        ch_versions = ch_versions.mix(MASK.out.versions)
        ch_amplicons = MASK.out.fasta
    }

    // Build reusable accession-to-taxid lookups for each database.
    BUILD_DB_TAXIDS(ch_dbs, ch_accession_taxonomy)
    ch_versions = ch_versions.mix(BUILD_DB_TAXIDS.out.versions)

    // Cluster amplicons and parse cluster membership rows.
    VSEARCH_CLUSTER_FAST(ch_amplicons)
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_FAST.out.versions)

    PARSE_UC(VSEARCH_CLUSTER_FAST.out.uc)
    ch_versions = ch_versions.mix(PARSE_UC.out.versions)

    // Add taxids to clustered accessions.
    JOIN_ACCESSION_TAXONOMY(
        combine_by_db(PARSE_UC.out.tsv, BUILD_DB_TAXIDS.out.accession_taxid)
    )
    ch_versions = ch_versions.mix(JOIN_ACCESSION_TAXONOMY.out.versions)

    CLUSTER_CONSENSUS(
        JOIN_ACCESSION_TAXONOMY.out.tsv,
        ch_taxdump,
    )
    ch_versions = ch_versions.mix(CLUSTER_CONSENSUS.out.versions)

    CLUSTER_CONSENSUS_REPORTING(CLUSTER_CONSENSUS.out.tsv)
    ch_versions = ch_versions.mix(CLUSTER_CONSENSUS_REPORTING.out.versions)
    multiqc_files = multiqc_files.mix(
        CLUSTER_CONSENSUS_REPORTING.out.resolution,
        CLUSTER_CONSENSUS_REPORTING.out.summary,
        CLUSTER_CONSENSUS_REPORTING.out.cluster_size,
        CLUSTER_CONSENSUS_REPORTING.out.taxon_search,
        // CLUSTER_CONSENSUS_REPORTING.out.amplified_taxa, // Disabled: makes reports slow.
    )

    // Summarize taxid distribution for each database.
    DB_DISTRIBUTION(
        BUILD_DB_TAXIDS.out.taxids_counts,
        ch_taxdump,
    )
    ch_versions = ch_versions.mix(DB_DISTRIBUTION.out.versions)

    if (params.taxon) {
        TAXONOMIC_COVERAGE(
            combine_by_db(CLUSTER_CONSENSUS.out.tsv, BUILD_DB_TAXIDS.out.taxids),
            params.taxon,
        )
        ch_versions = ch_versions.mix(TAXONOMIC_COVERAGE.out.versions)
        multiqc_files = multiqc_files.mix(
            TAXONOMIC_COVERAGE.out.mqc_summary,
            TAXONOMIC_COVERAGE.out.mqc_details,
        )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Create the settings table and MultiQC reports.
    REPORTING(
        multiqc_files,
        ch_primers_with_db,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.mix(BUILD_DB_TAXIDS.out.mqc.map { _meta, f -> f }),
    )

    emit:
    qc = REPORTING.out.html
    consensus = CLUSTER_CONSENSUS.out.tsv
}
