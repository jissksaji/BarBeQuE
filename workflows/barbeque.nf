/*
Import modules
*/
include { INPUT_CHECK } from './../modules/input_check'
include { PRIMER_SET } from './../subworkflows/primer_set'
include { COLLAPSE_PRIMERS as COLLAPSE_INPUT_PRIMERS } from './../modules/helper/collapse_primers'
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
include { FILTER_ACCESSIONS } from './../modules/helper/filter_accessions'
include { HIERARCHICAL_CLUSTERING } from './../subworkflows/hierarchical_clustering'
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
    ch_multiqc_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : channel.value([])
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect() : channel.value([])
    ch_versions = ch_db_versions
    multiqc_files = channel.from([])
    // The pre-installed (or user-supplied) taxdump folder used for taxonomic lookups
    def taxdump_path = params.taxdump ?: params.references.taxdump
    ch_taxdump = channel.value(
        file(taxdump_path, checkIfExists: !params.build_references)
    )

    // nucl_gb.accession2taxid is installed beside new_taxdump in taxonomy/.
    def accession_file = params.accession_taxonomy
    if (!accession_file && params.reference_base) {
        def taxonomy_dir = "${params.reference_base}/barbeque/${params.reference_version}/taxonomy"
        accession_file = files("${taxonomy_dir}/*{accession2taxid,genbank2taxid,nucl_}*").find()
    }
    // Compatibility with installations made before the shared taxonomy folder.
    if (!accession_file && params.reference_base) {
        def genbank_dir = "${params.reference_base}/barbeque/${params.reference_version}/genbank2taxid"
        accession_file = files("${genbank_dir}/*{accession2taxid,genbank2taxid,nucl_}*").find()
    }
    if (!accession_file && taxdump_path) {
        accession_file = files("${taxdump_path}/*{accession2taxid,genbank2taxid,nucl_}*").find()
    }

    if (!accession_file) {
        log.error("No accession-to-taxonomy mapping found - provide --accession_taxonomy <file>, or rebuild references so taxonomy/ contains both new_taxdump and nucl_gb.accession2taxid")
        System.exit(1)
    }

    ch_accession_taxonomy = channel.value(
        file(accession_file, checkIfExists: true)
    )



    // Primers come from either a named --primer_set (resolved live from the FooDMe2 catalog),
    // a hand-written --input samplesheet, or an --input directory of primer FASTAs - mutually exclusive, enforced in
    // lib/WorkflowPipeline.groovy. All converge on the same [primer, fwd, rev, min, max]
    // meta shape, so nothing below this point needs to know which path was used.
    if (params.primer_set) {
        PRIMER_SET()
        ch_primers = PRIMER_SET.out.primers
        ch_versions = ch_versions.mix(PRIMER_SET.out.versions)
    }
    else {
        def input_path = file(params.input, checkIfExists: true)
        if (input_path.isDirectory()) {
            ch_primer_fasta = channel.fromPath("${params.input}/*.{fa,fasta,fna}", checkIfExists: true)
                .map { fasta ->
                    def primer_id = fasta.baseName.replaceAll(/[\s\/\\:*?"<>|]/, '_')
                    tuple(
                        [id: primer_id, min: params.primer_min, max: params.primer_max],
                        fasta,
                    )
                }

            // Collapse each FASTA's fwd/rev-labelled variants into one consensus fwd + rev pair.
            COLLAPSE_INPUT_PRIMERS(ch_primer_fasta)
            ch_versions = ch_versions.mix(COLLAPSE_INPUT_PRIMERS.out.versions)

            ch_primers = COLLAPSE_INPUT_PRIMERS.out.fasta.map { meta, fasta ->
                def seqs = [:]
                def id = null
                def seq = new StringBuilder()
                fasta.eachLine { line ->
                    if (line.startsWith('>')) {
                        if (id) {
                            seqs[id] = seq.toString()
                        }
                        id = line.substring(1).trim()
                        seq = new StringBuilder()
                    }
                    else {
                        seq.append(line.trim())
                    }
                }
                if (id) {
                    seqs[id] = seq.toString()
                }
                [
                    primer: meta.id,
                    fwd: seqs["${meta.id}_fwd"],
                    rev: seqs["${meta.id}_rev"],
                    min: meta.min,
                    max: meta.max,
                ]
            }
        }
        else {
            samplesheet = channel.fromPath(input_path)

            // Check if the samplesheet is valid
            INPUT_CHECK(samplesheet)

            // Copy the samplesheet to the results folder
            STAGE_SAMPLESHEET(samplesheet)

            ch_primers = INPUT_CHECK.out.primers
        }
    }

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
    if (!(insilico_tool in ['obipcr', 'cutadapt'])) {
        log.error("Invalid --insilico_tool '${insilico_tool}' - must be 'obipcr' or 'cutadapt'")
        System.exit(1)
    }

    if (insilico_tool == 'cutadapt') {
        // Estimate a read buffer size from the largest sequence in each database
        COMPUTE_BUFFER(ch_dbs)
        // Run in-silico PCR with cutadapt to extract amplicons for each primer/database pair
        CUTADAPT_INSILICOPCR(
            ch_primers_with_db.map { meta, db -> tuple(meta.db, meta, db) }.combine(
                COMPUTE_BUFFER.out.buffersize.map { meta, buffersize -> tuple(meta.id, buffersize) },
                by: 0
            ).map { _db, meta, db, buffersize -> tuple(meta, db, buffersize) }
        )
        ch_versions = ch_versions.mix(COMPUTE_BUFFER.out.versions, CUTADAPT_INSILICOPCR.out.versions)
        ch_insilico_fasta = CUTADAPT_INSILICOPCR.out.fasta
    }
    else {
        // Run in-silico PCR with obipcr to extract amplicons for each primer/database pair
        OBIPCR_INSILICOPCR(
            ch_primers_with_db
        )
        ch_versions = ch_versions.mix(OBIPCR_INSILICOPCR.out.versions)
        ch_insilico_fasta = OBIPCR_INSILICOPCR.out.fasta

        // Reformat obipcr's raw output into a standard amplicon FASTA
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
        // Mask primer-binding regions in each amplicon to mimic real sequencing reads
        MASK(ch_insilico_by_status.valid)
        ch_amplicons = MASK.out.fasta
        ch_versions = ch_versions.mix(MASK.out.versions)
    }
    else {
        ch_amplicons = ch_insilico_by_status.valid
    }

    // Record the amplicon length distribution per primer/database (on the raw amplicons)
    AMPLICON_LENGTH(ch_amplicons)
    ch_versions = ch_versions.mix(AMPLICON_LENGTH.out.versions)
    multiqc_files = multiqc_files.mix(AMPLICON_LENGTH.out.tsv)

    // Map each database's accessions to taxids once per db (expensive full scan of genbank2taxid).
    // Moved ahead of clustering so amplicon accessions can be resolved to species for the divergence
    // screen below; still runs on the whole DB, so all downstream consumers are unchanged.
    BUILD_DB_TAXIDS(ch_dbs, ch_accession_taxonomy)
    ch_versions = ch_versions.mix(BUILD_DB_TAXIDS.out.versions)

    // Advisory per-species divergence screening (hierarchical clustering) - flags candidate
    // misannotations for review but does not change what gets clustered (see
    // params.exclude_accessions for the actual filter).
    if (params.screen_species_divergence) {
        HIERARCHICAL_CLUSTERING(ch_amplicons, BUILD_DB_TAXIDS.out.accession_taxid, ch_taxdump)
        ch_versions = ch_versions.mix(HIERARCHICAL_CLUSTERING.out.versions)
    }

    // Optionally drop author-curated bad accessions before clustering (the only step that changes
    // clustering input); otherwise cluster the amplicons unchanged.
    if (params.exclude_accessions) {
        FILTER_ACCESSIONS(ch_amplicons, file(params.exclude_accessions, checkIfExists: true))
        ch_to_cluster = FILTER_ACCESSIONS.out.fasta
        ch_versions = ch_versions.mix(FILTER_ACCESSIONS.out.versions)
    }
    else {
        ch_to_cluster = ch_amplicons
    }

    // Cluster amplicons into OTUs
    VSEARCH_CLUSTER_FAST(ch_to_cluster)
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER_FAST.out.versions)

    // Extract accession-to-cluster assignments from the vsearch .uc file
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

    if (params.completeness_table) {
        // Report how completely the target taxon is represented in each database
        COMPLETENESS_TABLE(
            BUILD_DB_TAXIDS.out.taxids_counts
        )
        ch_versions = ch_versions.mix(COMPLETENESS_TABLE.out.versions)
    }

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
