/*
Import modules
*/
include { INPUT_CHECK } from './../modules/input_check'
include { PRIMER_SET } from './../subworkflows/primer_set'
include { PARSE_PRIMERS } from './../modules/parse_primers/main'
include { MULTIQC } from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'
include { STAGE_FILE as STAGE_SAMPLESHEET } from './../modules/helper/stage_file'

/*
Helper Modules
*/
include { BUILD_DB_TAXIDS } from './../modules/helper/build_db_taxids'
include { DB_DISTRIBUTION } from './../modules/db_distribution/main.nf'
include { AMPLICON_LENGTH } from './../modules/seqkit/amplicon_lengths'
include { SPECIES_REPRESENTATION } from './../modules/helper/species_representation/main'
include { MASK } from './../modules/mask/main.nf'
include { PARSE_UC } from './../modules/helper/parse_uc'
include { JOIN_ACCESSION_TAXONOMY } from './../modules/helper/join_accession_taxonomy'
include { ACCESSION_BLOCKLIST } from './../modules/helper/accession_blocklist'
/*
Core Modules
*/
include { OBIPCR_INSILICOPCR } from './../modules/obipcr'
include { PARSE_OBIPCR } from './../modules/parse_obipcr/main.nf'
include { VSEARCH_CLUSTER_FAST } from './../modules/vsearch/cluster_fast'
include { VSEARCH_DEREPLICATION } from './../modules/vsearch/dereplication'
include { TAXONOMIC_COVERAGE } from './../modules/helper/taxonomic_coverage'
include { CLUSTER_CONSENSUS } from './../modules/helper/cluster_consensus'



workflow BARBEQUE {
    take:
    ch_dbs
    ch_db_versions

    main:
    ch_multiqc_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : channel.value([])
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : channel.value([])
    ch_versions = ch_db_versions
    multiqc_files = channel.from([])
    // The pre-installed (or user-supplied) taxdump folder used for taxonomic lookups
    def taxdump_path = params.taxdump ?: params.references.taxdump
    ch_taxdump = channel.value(
        file(taxdump_path, checkIfExists: true)
    )

    // nucl_gb.accession2taxid is installed beside new_taxdump in taxonomy/.
    def accession_file = params.accession_taxonomy
    if (!accession_file && params.reference_base) {
        def taxonomy_dir = "${params.reference_base}/barbeque/${params.reference_version}/taxonomy"
        accession_file = files("${taxonomy_dir}/*{accession2taxid,genbank2taxid,nucl_}*").find()
    }

    if (!accession_file) {
        log.error("No accession-to-taxonomy mapping found - provide --accession_taxonomy <file>, or rebuild references so taxonomy/ contains both new_taxdump and nucl_gb.accession2taxid")
        System.exit(1)
    }

    ch_accession_taxonomy = channel.value(
        file(accession_file, checkIfExists: true)
    )



    // Primers come from either a named --primer_set (resolved live from the FooDMe2 catalog),
    // a hand-written --input samplesheet, a single --input primer FASTA, or an --input directory
    // of primer FASTAs - mutually exclusive, enforced in lib/WorkflowPipeline.groovy. All converge
    // on the same [primer, fwd, rev, min, max] meta shape, so nothing below this point needs to
    // know which path was used.
    if (params.primer_set) {
        PRIMER_SET()
        samplesheet = PRIMER_SET.out.samplesheet
        ch_versions = ch_versions.mix(PRIMER_SET.out.versions)
    }
    else if (WorkflowPipeline.isFastaInput(params.input)) {
        // FASTA input (one file or a folder of them) is converted into the same samplesheet a
        // user would write by hand, so every route shares the validation and staging below.
        PARSE_PRIMERS(
            channel.value([
                [id: 'primers', min: params.primer_min, max: params.primer_max],
                file(params.input, checkIfExists: true),
            ])
        )
        ch_versions = ch_versions.mix(PARSE_PRIMERS.out.versions)
        samplesheet = PARSE_PRIMERS.out.samplesheet.map { _meta, tsv -> tsv }
    }
    else {
        samplesheet = channel.fromPath(file(params.input, checkIfExists: true))
    }

    // Check if the samplesheet is valid
    INPUT_CHECK(samplesheet)

    // Copy the samplesheet to the results folder
    STAGE_SAMPLESHEET(samplesheet)

    ch_primers = INPUT_CHECK.out.primers

    /*
     Combine each primer set with all requested databases
     [ meta, database_meta, database_path ]
    */
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



    def insilico_tool = (params.insilico_tool ?: 'obipcr').toLowerCase()
    if (insilico_tool != 'obipcr') {
        log.error("Invalid --insilico_tool '${insilico_tool}' - must be 'obipcr'")
        System.exit(1)
    }

    // Run in-silico PCR with obipcr to extract amplicons for each primer/database pair
    OBIPCR_INSILICOPCR(
        ch_primers_with_db
    )
    ch_versions = ch_versions.mix(OBIPCR_INSILICOPCR.out.versions)

    // Reformat obipcr's raw output into a standard amplicon FASTA
    PARSE_OBIPCR(
        OBIPCR_INSILICOPCR.out.raw_fasta
    )
    ch_versions = ch_versions.mix(PARSE_OBIPCR.out.versions)

    // Optionally remove listed accessions after OBI-PCR has been parsed. The filtered
    // FASTA is the only amplicon channel exposed to masking, length summaries, clustering,
    // taxonomy, and reports; the paired filtered TSV remains the published parsed result.
    if (params.accession_blocklist) {
        // The two processes emit independently. Attach the same scalar pair key
        // to both channels so join cannot depend on file arrival order or map identity.
        ch_obipcr_fasta_by_key = OBIPCR_INSILICOPCR.out.fasta.map { meta, fasta ->
            tuple("${meta.primer}|${meta.db}", meta, fasta)
        }
        ch_obipcr_tsv_by_key = PARSE_OBIPCR.out.tsv.map { meta, tsv ->
            tuple("${meta.primer}|${meta.db}", tsv)
        }
        ch_obipcr_with_tsv = ch_obipcr_fasta_by_key
            .join(ch_obipcr_tsv_by_key)
            // Drop the temporary join key; downstream modules use the original meta map.
            .map { _key, meta, fasta, tsv -> tuple(meta, fasta, tsv) }

        ACCESSION_BLOCKLIST(
            ch_obipcr_with_tsv,
            channel.value(file(params.accession_blocklist, checkIfExists: true)),
        )
        ch_versions = ch_versions.mix(ACCESSION_BLOCKLIST.out.versions)
        ch_insilico_fasta = ACCESSION_BLOCKLIST.out.fasta
    }
    else {
        ch_insilico_fasta = OBIPCR_INSILICOPCR.out.fasta
    }

    ch_insilico_fasta
        .filter { _meta, fasta -> file(fasta).size() > 0 }
        .set { ch_nonempty_insilico_fasta }

    if (params.dereplicate_amplicons) {
        VSEARCH_DEREPLICATION(ch_nonempty_insilico_fasta)
        ch_versions = ch_versions.mix(VSEARCH_DEREPLICATION.out.versions)
        ch_pre_mask_amplicons = VSEARCH_DEREPLICATION.out.fasta
    }
    else {
        ch_pre_mask_amplicons = ch_nonempty_insilico_fasta
    }

    if (params.mask) {
        // Mask primer-binding regions in each amplicon to mimic real sequencing reads
        MASK(ch_pre_mask_amplicons)
        ch_amplicons = MASK.out.fasta
        ch_versions = ch_versions.mix(MASK.out.versions)
    }
    else {
        ch_amplicons = ch_pre_mask_amplicons
    }

    AMPLICON_LENGTH(ch_amplicons)
    ch_versions = ch_versions.mix(AMPLICON_LENGTH.out.versions)
    multiqc_files = multiqc_files.mix(AMPLICON_LENGTH.out.tsv)

    // Map each database's accessions to taxids once per db (expensive full scan of genbank2taxid).
    BUILD_DB_TAXIDS(ch_dbs, ch_accession_taxonomy)
    ch_versions = ch_versions.mix(BUILD_DB_TAXIDS.out.versions)

    // Cluster the retained amplicons into OTUs.
    VSEARCH_CLUSTER_FAST(ch_amplicons)
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_FAST.out.versions)

    PARSE_UC(VSEARCH_CLUSTER_FAST.out.uc)
    ch_versions = ch_versions.mix(PARSE_UC.out.versions)

    // Join each cluster's accessions to their taxids using the per-db lookup
    JOIN_ACCESSION_TAXONOMY(
        PARSE_UC.out.tsv.map { m, tsv -> tuple(m.db, m, tsv) }.combine(
            BUILD_DB_TAXIDS.out.accession_taxid.map { m, f -> tuple(m.id, f) },
            by: 0
        ).map { _db_id, m, tsv, f -> tuple(m, tsv, f) }
    )
    ch_versions = ch_versions.mix(JOIN_ACCESSION_TAXONOMY.out.versions)

    // Derive a consensus taxonomic assignment for each cluster
    CLUSTER_CONSENSUS(
        JOIN_ACCESSION_TAXONOMY.out.tsv,
        ch_taxdump,
    )
    ch_versions = ch_versions.mix(CLUSTER_CONSENSUS.out.versions)
    multiqc_files = multiqc_files.mix(CLUSTER_CONSENSUS.out.tsv)

    // Summarise the taxonomic distribution of each reference database
    DB_DISTRIBUTION(
        BUILD_DB_TAXIDS.out.taxids_counts,
        ch_taxdump,
    )
    ch_versions = ch_versions.mix(DB_DISTRIBUTION.out.versions)

    if (params.taxon) {
        // Assess taxonomic coverage of clusters against the requested taxon
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
        multiqc_files = multiqc_files.mix(TAXONOMIC_COVERAGE.out.tsv)

        // Combine coverage and consensus results into a per-species representation table
        SPECIES_REPRESENTATION(
            TAXONOMIC_COVERAGE.out.tsv.join(CLUSTER_CONSENSUS.out.tsv, by: 0)
        )
        ch_versions = ch_versions.mix(SPECIES_REPRESENTATION.out.versions)
        multiqc_files = multiqc_files.mix(SPECIES_REPRESENTATION.out.tsv)
    }

    // Collate software versions from all executed processes
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Combine by meta dict to generate separate reports for each primer-db combination
    multiqc_by_set = multiqc_files.groupTuple(by: 0)

    // Generate one MultiQC report per primer-database combination
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
