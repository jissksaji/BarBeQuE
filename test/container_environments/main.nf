nextflow.enable.dsl = 2

process CHECK_ENVIRONMENT {
    tag "${module_name}"

    conda { environment_file }

    input:
    tuple val(module_name), val(environment_file), val(check_command)

    output:
    tuple val(module_name), path("${module_name}.ok"), emit: results

    script:
    """
    set -euo pipefail
    ${check_command}
    touch ${module_name}.ok
    """
}

workflow {
    def module_root = file("${projectDir}/../../modules", checkIfExists: true)
    def checks = [
        ['custom_dumpsoftwareversions', "${module_root}/custom/dumpsoftwareversions/environment.yml", 'multiqc --version && python -c "import yaml"'],
        ['db_distribution', "${module_root}/db_distribution/environment.yml", 'python -c "import numpy, pandas, plotly, taxidTools" && sed --version'],
        ['download', "${module_root}/download/environment.yml", 'bash --version && wget --version && gzip --version && sha256sum --version && grep --version && python --version'],
        ['gunzip', "${module_root}/gunzip/environment.yml", 'gunzip --version && sed --version'],
        ['accession_blocklist', "${module_root}/helper/accession_blocklist/environment.yml", 'python3 --version && sed --version'],
        ['build_db_taxids', "${module_root}/helper/build_db_taxids/environment.yml", 'awk --version && grep --version && sed --version && sort --version'],
        ['cluster_consensus', "${module_root}/helper/cluster_consensus/environment.yml", 'python3 -c "import taxidTools" && sed --version'],
        ['join_accession_taxonomy', "${module_root}/helper/join_accession_taxonomy/environment.yml", 'awk --version && head --version'],
        ['parse_uc', "${module_root}/helper/parse_uc/environment.yml", 'awk --version && head --version'],
        ['species_representation', "${module_root}/helper/species_representation/environment.yml", 'python3 --version && sed --version'],
        ['stage_file', "${module_root}/helper/stage_file/environment.yml", 'touch --version'],
        ['taxid_db_filter', "${module_root}/helper/taxid_db_filter/environment.yml", 'taxonkit version && python3 -c "import Bio" && awk --version && grep --version && sed --version'],
        ['taxonomic_coverage', "${module_root}/helper/taxonomic_coverage/environment.yml", 'python3 -c "import ete3, pandas, plotly" && sed --version'],
        ['mask', "${module_root}/mask/environment.yml", 'python3 --version && sed --version'],
        ['multiqc', "${module_root}/multiqc/environment.yml", 'multiqc --version && sed --version'],
        ['obipcr', "${module_root}/obipcr/environment.yml", 'obipcr --version && seqkit version && sed --version && head --version'],
        ['parse_obipcr', "${module_root}/parse_obipcr/environment.yml", 'python --version && sed --version'],
        ['parse_primers', "${module_root}/parse_primers/environment.yml", 'python3 --version && sed --version'],
        ['amplicon_lengths', "${module_root}/seqkit/amplicon_lengths/environment.yml", 'seqkit version && sed --version'],
        ['db_filter', "${module_root}/seqkit/db_filter/environment.yml", 'seqkit version && sed --version'],
        ['streamlit', "${module_root}/streamlit/environment.yml", 'python -m streamlit --version && python -c "import matplotlib, pandas, plotly, venn" && pkill --version && sed --version'],
        ['untar', "${module_root}/untar/environment.yml", 'tar --version && grep --version && sed --version && rm --version'],
        ['vsearch_cluster_fast', "${module_root}/vsearch/cluster_fast/environment.yml", 'vsearch --version && sed --version && head --version'],
        ['vsearch_dereplication', "${module_root}/vsearch/dereplication/environment.yml", 'vsearch --version && sed --version && head --version'],
    ]

    CHECK_ENVIRONMENT(channel.fromList(checks))
    CHECK_ENVIRONMENT.out.results.view { module_name, _marker -> "CONTAINER_OK=${module_name}" }
}
