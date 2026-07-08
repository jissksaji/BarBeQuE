/*
Include Modules
*/

include { DOWNLOAD_AND_COMBINE_DB } from './../modules/download/main'
include { UNTAR as UNTAR_TAXDUMP } from './../modules/untar'
include { GUNZIP as GUNZIP_GENBANK } from './../modules/gunzip'

workflow BUILD_REFERENCES {

    def dbs_to_download = [
        [[db: "plastid", dir: "refseq/plastid", filename: "plastid.genomic.fna"], ["https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.1.1.genomic.fna.gz", "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.2.1.genomic.fna.gz", "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.3.1.genomic.fna.gz"]],
        [[db: "mitochondrion", dir: "refseq/mitochondrion", filename: "mitochondrion.genomic.fna"], ["https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz"]],
        [[db: "plasmid", dir: "refseq/plasmid", filename: "plasmid.genomic.fna"], [
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.1.1.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.1.2.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.2.1.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.2.2.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.3.1.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.3.2.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.4.1.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.4.2.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.5.1.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.5.2.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.6.1.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.6.2.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.7.1.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.7.2.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.8.1.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.8.2.genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.9.1.genomic.fna.gz",
        ]],
    ]

    // Every database in conf/resources.config that declares a `source` is fetched straight
    // from params.reference_sources (conf/reference_sources.config), publishing to the same
    // path already declared as its `db` value - so no directory/filename layout is duplicated.
    def path_prefix = "${params.reference_base}/barbeque/${params.reference_version}/"
    params.references.databases.each { db_id, info ->
        if (info.source) {
            def source = params.reference_sources[info.source]
            def target = file(info.db)
            dbs_to_download << [
                [db: db_id, dir: (target.parent.toString() - path_prefix), filename: target.name],
                [source.url],
            ]
        }
    }

    DOWNLOAD_AND_COMBINE_DB(channel.fromList(dbs_to_download))

    // Both large, optional, off by default - see --install_taxdump / --install_genbank
    if (params.install_taxdump) {
        UNTAR_TAXDUMP(
            channel.of([[id: 'new_taxdump'], file(params.reference_sources[params.references.taxdump_source].url)])
        )
    }

    if (params.install_genbank) {
        GUNZIP_GENBANK(
            channel.of([[id: 'nucl_gb.accession2taxid'], file(params.reference_sources.ncbi_genbank_accession2taxid.url)])
        )
    }

    if (params.build_references) {
        workflow.onComplete = {
            log.info('Installation complete - deleting staged files. ')
            log.info("Reference files at ${params.reference_base}")
            workDir.resolve("stage-${workflow.sessionId}").deleteDir()
        }
    }
}
