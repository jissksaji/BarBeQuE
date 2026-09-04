process STAGE_FILE {
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    path(f)

    output:
    path(f)                 , emit: staged_file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    touch dummy.txt $args
    
    """
}
