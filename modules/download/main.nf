process DOWNLOAD_AND_COMBINE_DB {
    tag "${meta.db}"
    label 'medium_serial'

    input:
    tuple val(meta), val(urls)

    output:
    tuple val(meta), path("${meta.filename}"), emit: db, optional: true
    path 'versions.yml', emit: versions

    script:
    def urlList = urls.collect { "\"${it}\"" }.join(' ')

    """
    mkdir downloads

    # Download every part of this database. Unavailable parts are retried by
    # wget and then skipped so that one broken upstream link does not stop the
    # remaining reference downloads.
    for url in ${urlList}; do
        name=\$(basename "\$url")
        tmp="downloads/.\${name}.part"
        if [[ "\$url" == file://* ]]; then
            if ! cp "\${url#file://}" "\$tmp"; then
                echo "WARNING: Could not copy \$url; skipping it." >&2
                rm -f "\$tmp"
                continue
            fi
        else
            downloaded=false
            for attempt in 1 2 3; do
                if wget --tries=1 --timeout=60 -O "\$tmp" "\$url"; then
                    downloaded=true
                    break
                fi
                echo "WARNING: Download attempt \$attempt failed for \$url." >&2
                rm -f "\$tmp"
            done
            if [[ "\$downloaded" != true ]]; then
                echo "WARNING: Could not download \$url after retries; skipping it." >&2
                continue
            fi
        fi
        mv "\$tmp" "downloads/\$name"
    done

    # Decompress and combine all parts into one file.
    shopt -s nullglob
    for part in downloads/*; do
        if [[ "\$part" == *.gz ]]; then
            gzip -cd "\$part" >> "${meta.filename}"
        else
            cat "\$part" >> "${meta.filename}"
        fi
    done

    # MetaFish is supplied as CSV, so convert it to FASTA.
    if [[ -s "${meta.filename}" && "${meta.format}" == "metafish_csv" ]]; then
        mv "${meta.filename}" metafish.csv
        metafish_to_fasta.py metafish.csv "${meta.filename}"
    fi

    # An optional output lets Nextflow omit this database when every link
    # failed, while continuing to download all other databases.
    if [[ ! -s "${meta.filename}" ]] || ! grep -q '^>' "${meta.filename}"; then
        echo "WARNING: No valid FASTA was produced for ${meta.db}; skipping this database." >&2
        rm -f "${meta.filename}"
    fi

    if [[ -f "${meta.filename}" ]]; then
        output_sha256=\$(sha256sum "${meta.filename}" | cut -d' ' -f1)
    else
        output_sha256=skipped
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        source_release: '${meta.release}'
        output_sha256: \$output_sha256
        wget: \$(wget --version | head -n 1 | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    """
    printf '>stub\\nACGT\\n' > "${meta.filename}"
    echo '"${task.process}": {output_sha256: stub}' > versions.yml
    """
}
