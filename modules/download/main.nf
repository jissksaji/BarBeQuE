process DOWNLOAD_AND_COMBINE_DB {
    tag "${meta.db}"
    label 'medium_serial'

    input:
    tuple val(meta), val(urls)

    output:
    tuple val(meta), path("${meta.filename}"), emit: db
    path "versions.yml", emit: versions

    script:
    def url_list = urls instanceof List ? urls.join(' ') : urls
    """
    # Download all URLs
    for url in ${url_list}; do
        filename=\$(basename "\$url")
        wget -qO "\$filename" "\$url"
    done

    # Decompress and combine into a single file
    # We use zcat for .gz files and cat for uncompressed files
    for f in *; do
        if [[ "\$f" == *.gz ]]; then
            gunzip -c "\$f" >> "${meta.filename}"
        elif [[ "\$f" != "${meta.filename}" ]] && [[ "\$f" != "versions.yml" ]]; then
            cat "\$f" >> "${meta.filename}"
        fi
    done

    # Cleanup downloaded raw parts
    for f in *; do
        if [[ "\$f" != "${meta.filename}" ]] && [[ "\$f" != "versions.yml" ]]; then
            rm "\$f"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n 1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
