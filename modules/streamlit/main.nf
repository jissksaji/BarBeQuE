process STREAMLIT {

    tag "dashboard"
    label 'process_low'
    debug true
    cache false

    conda "${moduleDir}/environment.yml"

    input:
    val results_dir

    output:
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

pkill -f "streamlit run.*app.py" || true

python -m streamlit run ${moduleDir}/../../bin/app.py \\
    --server.port 8501 \\
    --server.headless true \\
    -- "${results_dir}" > "${params.outdir}/streamlit.log" 2>&1 &

sleep 3

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version | sed 's/Python //')
    streamlit: \$(python -m streamlit --version | sed 's/Streamlit, version //')
END_VERSIONS
    """
}
