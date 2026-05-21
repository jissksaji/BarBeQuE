params.outdir = "../results"

include { STREAMLIT } from '../modules/streamlit'

workflow INTERACTIVE_RESULTS {
    take:
    consensus_dir

    main:
    consensus_dir.map { it ->
        log.info("\n\033[1;32mPort 8501 is automatically forwarded. View the dashboard at http://localhost:8501 \nTo kill the application, use this terminal (Ctrl+C).\nStreamlit is unstable,might have to wait till all the processes are completed to visualise the data\033[0m\n")
        return it
    } | STREAMLIT
    STREAMLIT.out.versions.subscribe { }

    emit:
    versions = STREAMLIT.out.versions
}


workflow {
    INTERACTIVE_RESULTS(file("${params.outdir}/consensus"))
}
