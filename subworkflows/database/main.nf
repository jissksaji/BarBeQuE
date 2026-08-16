/*
Include Modules
*/
include { CUSTOM_DB_FILTER } from './../../modules/seqkit/custom_db_filter'
include { TAXID_DB_FILTER } from './../../modules/helper/taxid_db_filter'
include { BLOCKLIST_FILTER } from './../../modules/helper/blocklist_filter'

def has_fasta_header(db) {
    def first_record = null
    db.withReader { reader ->
        first_record = reader.iterator().find { line -> line.trim() }
    }
    return first_record?.startsWith('>')
}

workflow DATABASE {

    main:
    ch_versions = channel.from([])

    // the database to use - either pre-installed or user-provided
    // Pre-installed can be a list, coma-separated:  db1,db2,db3
    these_dbs = []
    if (params.custom_db) {
        def custom_file = file(params.custom_db, checkIfExists: true)
        these_dbs << [["id": custom_file.baseName], custom_file]
    }
    else if (params.dbs) {
        valid_databases = params.references.databases.keySet()
        params.dbs
            .split(",")
            .collect { db_str -> db_str.toLowerCase() }
            .each { db ->
                if (!valid_databases.contains(db)) {
                    log.info("Not a valid database: ${db}\nValid options are: ${valid_databases}\n")
                    System.exit(1)
                }
                def info = params.references.databases[db]
                if (info.prebuilt) {
                    log.info("'${db}' is a pre-built BLAST index, not an in-silico-PCR reference - it cannot be used with --dbs.")
                    System.exit(1)
                }
                these_dbs << [
                    ["id": db],
                    file(info.db, checkIfExists: true),
                ]
            }
    }
    ch_dbs = channel.fromList(these_dbs)

    // Both taxonomy filters need an accession-to-taxid mapping. Resolve it once
    // here because DATABASE runs before the main BARBEQUE workflow.
    if (params.taxid || params.blocklist) {
        def accession_file = params.accession_taxonomy
        if (!accession_file && params.reference_base) {
            def taxonomy_dir = "${params.reference_base}/barbeque/${params.reference_version}/taxonomy"
            accession_file = files("${taxonomy_dir}/*{accession2taxid,genbank2taxid,nucl_}*").find()
        }
        if (!accession_file && params.reference_base) {
            def genbank_dir = "${params.reference_base}/barbeque/${params.reference_version}/genbank2taxid"
            accession_file = files("${genbank_dir}/*{accession2taxid,genbank2taxid,nucl_}*").find()
        }
        def taxdump_path = params.taxdump ?: params.references.taxdump
        if (!accession_file && taxdump_path) {
            accession_file = files("${taxdump_path}/*{accession2taxid,genbank2taxid,nucl_}*").find()
        }
        if (!accession_file) {
            log.error("--taxid and --blocklist require an accession-to-taxid mapping - provide --accession_taxonomy <file>, or rebuild the taxonomy references")
            System.exit(1)
        }
        ch_accession_taxonomy = channel.value(
            file(accession_file, checkIfExists: true)
        )
    }

    // Restrict every selected/custom db to a single taxon (and its descendants) before anything
    // else runs against it - lets --taxid be combined with --dbs or --custom_db to analyse/check
    // just that taxon instead of the whole database.
    if (params.taxid) {
        // Same taxdump/accession_taxonomy resolution as workflows/barbeque.nf, duplicated here
        // since DATABASE() runs before BARBEQUE() in main.nf and needs it to resolve --taxid.
        def taxdump_path = params.taxdump ?: params.references.taxdump
        ch_taxdump_for_taxid = channel.value(
            file(taxdump_path, checkIfExists: !params.build_references)
        )

        TAXID_DB_FILTER(ch_dbs, ch_taxdump_for_taxid, ch_accession_taxonomy, params.taxid)
        ch_versions = ch_versions.mix(TAXID_DB_FILTER.out.versions)
        ch_dbs = TAXID_DB_FILTER.out.fasta
    }

    // Optional FooDMe2 blocklist. The source is pinned in BlocklistCatalog.groovy,
    // just like named primer sources are pinned in PrimerCatalog.groovy.
    if (params.blocklist) {
        ch_blocklist = channel.value(
            file(BlocklistCatalog.url(), checkIfExists: true)
        )
        BLOCKLIST_FILTER(ch_dbs, ch_accession_taxonomy, ch_blocklist)
        ch_versions = ch_versions.mix(BLOCKLIST_FILTER.out.versions)
        ch_dbs = BLOCKLIST_FILTER.out.fasta
    }

    // filter the db, re-using a cached filtered copy from a previous run if present
    if (params.custom_db_filter) {
        ch_dbs
            .branch { meta, db ->
                // The cache is keyed only by meta.id, so taxonomy-filtered runs must never
                // read (or clobber) the cache of the whole, unfiltered database.
                def cached_db = params.reference_base ? file("${params.reference_base}/barbeque/${params.reference_version}/filtered/${meta.id}/${meta.id}.cleaned.fasta") : null
                cached: !params.taxid && !params.blocklist && cached_db?.exists() && has_fasta_header(cached_db)
                    return tuple(meta, cached_db)
                uncached: true
                    return tuple(meta, db)
            }
            .set { ch_dbs_by_cache }

        CUSTOM_DB_FILTER(ch_dbs_by_cache.uncached)
        ch_versions = ch_versions.mix(CUSTOM_DB_FILTER.out.versions)

        ch_dbs = ch_dbs_by_cache.cached.mix(CUSTOM_DB_FILTER.out.fasta)
    }

    emit:
    db       = ch_dbs
    versions = ch_versions
}
