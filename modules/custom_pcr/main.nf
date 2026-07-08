process CUSTOM_PCR_CPP {

    tag "${meta.primer}|${meta.db}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path('*_insilico.fasta'), emit: fasta
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.primer}_${meta.db}"

    def mismatches = params.containsKey('custom_pcr_mismatches')
        ? params.custom_pcr_mismatches as int
        : params.cutadapt_mismatches as int

    """
    g++ -O3 -std=c++17 -pthread -o pcr "${moduleDir}/pcr.cpp"

    ./pcr \\
        -i "${db}" \\
        -f "${meta.fwd}" \\
        -r "${meta.rev}" \\
        -m ${mismatches} \\
        --min-length ${meta.min} \\
        --max-length ${meta.max} \\
        --threads ${task.cpus} \\
        -o "${prefix}_insilico.fasta" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
        pcr.cpp: \$(g++ --version | head -n1)
END_VERSIONS
    """
}
