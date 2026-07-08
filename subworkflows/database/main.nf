/*
Include Modules
*/
include { CUSTOM_DB_FILTER } from './../../modules/seqkit/custom_db_filter'

workflow DATABASE {

    main:
    ch_versions = channel.from([])

    // the database to use - either pre-installed or user-provided
    // Pre-installed can be a list, coma-separated:  db1,db2,db3
    these_dbs = []
    if (params.custom_db) {
        these_dbs << [["id": "custom"], file(params.custom_db, checkIfExists: true)]
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
                these_dbs << [
                    ["id": db],
                    file(params.references.databases[db].db, checkIfExists: true),
                ]
            }
    }
    ch_dbs = channel.fromList(these_dbs)

    // filter the db, re-using a cached filtered copy from a previous run if present
    if (params.custom_db_filter) {
        ch_dbs
            .branch { meta, db ->
                def cached_db = file("${params.reference_base}/barbeque/${params.reference_version}/filtered/${meta.id}/${meta.id}.cleaned.fasta")
                cached: cached_db.exists()
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
