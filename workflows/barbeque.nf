include { INPUT_CHECK } from './../modules/input_check'
include { PRIMER_SET } from './../subworkflows/primer_set/main'
include { PARSE_PRIMERS } from './../modules/parse_primers/main'
include { MULTIQC } from './../modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions/main'
include { STAGE_FILE as STAGE_SAMPLESHEET } from './../modules/helper/stage_file/main'
include { BUILD_DB_TAXIDS } from './../modules/helper/build_db_taxids/main'
include { DB_DISTRIBUTION } from './../modules/db_distribution/main'
include { AMPLICON_LENGTH } from './../modules/seqkit/amplicon_lengths/main'
include { SPECIES_REPRESENTATION } from './../modules/helper/species_representation/main'
include { MASK } from './../modules/mask/main'
include { PARSE_UC } from './../modules/helper/parse_uc/main'
include { JOIN_ACCESSION_TAXONOMY } from './../modules/helper/join_accession_taxonomy/main'
include { ACCESSION_BLOCKLIST } from './../modules/helper/accession_blocklist/main'
include { OBIPCR_INSILICOPCR } from './../modules/obipcr/main'
include { PARSE_OBIPCR } from './../modules/parse_obipcr/main'
include { VSEARCH_CLUSTER_FAST } from './../modules/vsearch/cluster_fast/main'
include { VSEARCH_DEREPLICATION } from './../modules/vsearch/dereplication/main'
include { TAXONOMIC_COVERAGE } from './../modules/helper/taxonomic_coverage/main'
include { CLUSTER_CONSENSUS } from './../modules/helper/cluster_consensus/main'

/*
 * Pair a per-primer channel with a per-database lookup table.
 * Amplicon results carry the database name in meta.db, the lookups from
 * BUILD_DB_TAXIDS carry it in meta.id, so the two cannot be joined directly.
 * Re-key both to the database name, combine, then emit [meta, file, lookup].
 */
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
    // MULTIQC declares config and logo as path inputs, which cannot be null, so an
    // unset one is passed as an empty list and stages no file.
    ch_multiqc_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : channel.value([])
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : channel.value([])
    ch_versions = ch_db_versions
    multiqc_files = channel.empty()

    // Three routes in: a named --primer_set pulled from the FooDMe2 catalog, an --input
    // primer FASTA (or a directory of them), or an --input samplesheet. --input and
    // --primer_set are mutually exclusive, enforced in lib/WorkflowPipeline.groovy. All
    // three end up as a samplesheet, so nothing below needs to know which was used.
    if (params.primer_set) {
        PRIMER_SET()
        samplesheet = PRIMER_SET.out.samplesheet
        ch_versions = ch_versions.mix(PRIMER_SET.out.versions)
    }
    else if (WorkflowPipeline.isFastaInput(params.input)) {
        // Turn the FASTA into the samplesheet a user would otherwise write by hand.
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

    // Requires every column, and sanitises primer names by replacing path-unsafe
    // characters with '_' - sets whose names collide once sanitised are rejected.
    INPUT_CHECK(samplesheet)

    // STAGE_FILE runs no real command; it re-emits its input so the publishDir rule in
    // conf/modules.config copies the samplesheet into pipeline_info/.
    STAGE_SAMPLESHEET(samplesheet)

    ch_primers = INPUT_CHECK.out.primers

    /*
     * One in-silico PCR per primer x database pair. combine() yields
     * [primer_meta, db_meta, db_fasta]; the map folds the database name into the primer
     * meta as .db, which is what every downstream tag and <primer>_<db> filename uses.
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


    OBIPCR_INSILICOPCR(
        ch_primers_with_db
    )
    ch_versions = ch_versions.mix(OBIPCR_INSILICOPCR.out.versions)

    // Parse the obipcr headers into a per-amplicon TSV. This has to read raw_fasta:
    // the _sub[start..stop] suffix and the trailing JSON that parse_obipcr.py needs are
    // both stripped from the .fasta output by seqkit inside the obipcr module.
    PARSE_OBIPCR(
        OBIPCR_INSILICOPCR.out.raw_fasta
    )
    ch_versions = ch_versions.mix(PARSE_OBIPCR.out.versions)

    // Optionally drop listed accessions. The filtered FASTA becomes the only amplicon
    // channel below, and conf/modules.config publishes the filtered TSV under
    // parsed_obipcr/ in place of the unfiltered one.
    if (params.accession_blocklist) {
        // The two processes emit independently, so pair them back up on the meta map
        // both carry unchanged from OBIPCR_INSILICOPCR; join matches it by value, not
        // by arrival order.
        ACCESSION_BLOCKLIST(
            OBIPCR_INSILICOPCR.out.fasta.join(PARSE_OBIPCR.out.tsv),
            channel.value(file(params.accession_blocklist, checkIfExists: true)),
        )
        ch_versions = ch_versions.mix(ACCESSION_BLOCKLIST.out.versions)
        ch_insilico_fasta = ACCESSION_BLOCKLIST.out.fasta
    }
    else {
        ch_insilico_fasta = OBIPCR_INSILICOPCR.out.fasta
    }

    // A primer that amplifies nothing in a database is a valid result, but an empty
    // FASTA has nothing for the steps below to work on, so drop the pair here.
    ch_insilico_fasta
        .filter { _meta, fasta -> fasta.size() > 0 }
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
        // Cut each amplicon down to what a real run would return: paired-end keeps
        // --read_length bases at each end and Ns out the unsequenced middle, single-end
        // truncates to --read_length. Amplicons the reads already span are left alone.
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

    // Resolve every accession in each database to a taxid. This is a full scan of the
    // accession-to-taxid table, so it runs once per database and the result is reused by
    // the per-primer joins below rather than being rebuilt for each primer.
    BUILD_DB_TAXIDS(ch_dbs, ch_accession_taxonomy)
    ch_versions = ch_versions.mix(BUILD_DB_TAXIDS.out.versions)

    // Only the .uc membership table is kept - the centroid FASTA output is disabled in
    // the module - so PARSE_UC turns it into cluster/accession rows.
    VSEARCH_CLUSTER_FAST(ch_amplicons)
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_FAST.out.versions)

    PARSE_UC(VSEARCH_CLUSTER_FAST.out.uc)
    ch_versions = ch_versions.mix(PARSE_UC.out.versions)

    // Join each cluster's accessions to their taxids using the per-db lookup
    JOIN_ACCESSION_TAXONOMY(
        combine_by_db(PARSE_UC.out.tsv, BUILD_DB_TAXIDS.out.accession_taxid)
    )
    ch_versions = ch_versions.mix(JOIN_ACCESSION_TAXONOMY.out.versions)

    CLUSTER_CONSENSUS(
        JOIN_ACCESSION_TAXONOMY.out.tsv,
        ch_taxdump,
    )
    ch_versions = ch_versions.mix(CLUSTER_CONSENSUS.out.versions)
    multiqc_files = multiqc_files.mix(CLUSTER_CONSENSUS.out.tsv)

    // Per-database summary, driven by the taxid counts rather than by any primer's run.
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
        multiqc_files = multiqc_files.mix(TAXONOMIC_COVERAGE.out.tsv)

        SPECIES_REPRESENTATION(
            TAXONOMIC_COVERAGE.out.tsv.join(CLUSTER_CONSENSUS.out.tsv, by: 0)
        )
        ch_versions = ch_versions.mix(SPECIES_REPRESENTATION.out.versions)
        multiqc_files = multiqc_files.mix(SPECIES_REPRESENTATION.out.tsv)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Group every collected table by its meta map, so MULTIQC emits one report
    // per primer-database combination rather than one for the whole run.
    MULTIQC(
        multiqc_files.groupTuple(by: 0),
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
        ch_multiqc_config,
        ch_multiqc_logo,
    )

    emit:
    qc = MULTIQC.out.html
    consensus = CLUSTER_CONSENSUS.out.tsv
}
