/*
Include Modules
*/
include { COLLAPSE_PRIMERS } from './../../modules/helper/collapse_primers'

workflow PRIMER_SET {

    main:
    ch_versions = channel.from([])

    // Resolve each requested --primer_set name against the live FooDMe2 catalog (PrimerCatalog,
    // lib/PrimerCatalog.groovy) - same catalog --list_primers prints, always current, nothing
    // hand-maintained here that could drift out of sync with upstream.
    def catalog = PrimerCatalog.fetchCatalog()
    def valid_names = catalog.keySet()

    def these_primers = []
    params.primer_set
        .split(',')
        .collect { name -> name.trim() }
        .each { name ->
            if (!valid_names.contains(name)) {
                log.info("Not a valid primer set: ${name}\nValid options are: ${valid_names}\n")
                System.exit(1)
            }
            def sub = PrimerCatalog.fetchSubConfig(catalog[name].config)
            if (!sub.fasta) {
                log.info("Primer set '${name}' has no fasta reference in its upstream config - cannot use it.\n")
                System.exit(1)
            }
            these_primers << [
                [id: name, min: sub.min, max: sub.max],
                file(PrimerCatalog.fastaUrl(sub.fasta), checkIfExists: true),
            ]
        }
    ch_primer_fasta = channel.fromList(these_primers)

    // Collapse each set's fwd/rev-labelled variants into one consensus fwd + rev pair
    COLLAPSE_PRIMERS(ch_primer_fasta)
    ch_versions = ch_versions.mix(COLLAPSE_PRIMERS.out.versions)

    // Turn the collapsed 2-record (fwd/rev) fasta back into the same meta-map shape
    // modules/input_check.nf already produces (primer/fwd/rev/min/max), so workflows/barbeque.nf
    // doesn't need to change anything downstream of it.
    ch_primers = COLLAPSE_PRIMERS.out.fasta.map { meta, fasta ->
        def seqs = [:]
        def id = null
        def seq = new StringBuilder()
        fasta.eachLine { line ->
            if (line.startsWith('>')) {
                if (id) {
                    seqs[id] = seq.toString()
                }
                id = line.substring(1).trim()
                seq = new StringBuilder()
            }
            else {
                seq.append(line.trim())
            }
        }
        if (id) {
            seqs[id] = seq.toString()
        }
        [
            primer: meta.id,
            fwd: seqs["${meta.id}_fwd"],
            rev: seqs["${meta.id}_rev"],
            min: meta.min,
            max: meta.max,
        ]
    }

    emit:
    primers  = ch_primers
    versions = ch_versions
}
