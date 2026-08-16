include { DOWNLOAD_AND_COMBINE_DB } from './../modules/download/main'
include { UNTAR as UNTAR_TAXDUMP } from './../modules/untar'
include { GUNZIP as GUNZIP_GENBANK } from './../modules/gunzip'

def enabled(value) { value.toString().toBoolean() }

workflow BUILD_REFERENCES {

    // Download every database that has URLs in conf/resources.config.
    def downloads = []
    params.references.databases.each { id, db ->
        if (db.urls) {
            downloads << [
                [
                    db: id,
                    dir: "databases/${id}",
                    filename: "${id}.fasta",
                    format: db.format,
                    release: db.release,
                ],
                db.urls,
            ]
        }
    }

    // Install the primer FASTAs used by FooDMe2.
    PrimerCatalog.fetchFastas().each { name, url ->
        downloads << [
            [
                db: "primer_${name}",
                dir: 'primers',
                filename: name,
                format: 'fasta',
                release: "FooDMe2 ${PrimerCatalog.REVISION}",
            ],
            [url],
        ]
    }

    DOWNLOAD_AND_COMBINE_DB(channel.fromList(downloads))

    // Install both NCBI taxonomy files together.
    if (enabled(params.install_taxdump)) {
        UNTAR_TAXDUMP(channel.of([
            [id: 'new_taxdump'],
            file(params.reference_sources.ncbi_taxdump.url),
        ]))
        GUNZIP_GENBANK(channel.of([
            [id: 'nucl_gb.accession2taxid'],
            file(params.reference_sources.ncbi_genbank_accession2taxid.url),
        ]))
    }

    workflow.onComplete = {
        if (workflow.success) {
            log.info("References installed in ${params.reference_base}")
        }
    }
}
