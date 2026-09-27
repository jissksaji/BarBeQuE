nextflow.enable.dsl = 2

include { DB_FILTER } from '../../modules/seqkit/db_filter/main.nf'

/*
 * DB_FILTER builds its seqkit pipeline from params.db_filter_*, so one run can only
 * exercise one combination. run_cases.sh drives this workflow once per combination
 * and passes the record IDs that should survive as --expect.
 *
 * Without --expect the workflow only reports what survived, so a plain
 * `nextflow run test/db_filter/main.nf` still works from run_all_tests.sh.
 */
params.expect = null

workflow {
    def meta = [id: 'test_db']
    def reference = file("${projectDir}/reference.fasta")

    DB_FILTER(tuple(meta, reference))

    DB_FILTER.out.fasta
        .map { _meta, fasta ->
            fasta.readLines()
                .findAll { line -> line.startsWith('>') }
                .collect { line -> line.substring(1).split(/\s/)[0] }
                .sort()
        }
        .subscribe { kept ->
            if (params.expect == null) {
                log.info("DB_FILTER kept ${kept.size()} record(s): ${kept.join(', ')}")
                return
            }

            def expected = params.expect.toString().tokenize(',').collect { it.trim() }.findAll { it }.sort()

            if (kept != expected) {
                def missing = expected - kept
                def extra = kept - expected
                error(
                    "DB_FILTER kept the wrong records.\n" +
                    "  expected: ${expected.join(', ')}\n" +
                    "  actual  : ${kept.join(', ')}\n" +
                    "  wrongly dropped: ${missing.join(', ') ?: '(none)'}\n" +
                    "  wrongly kept   : ${extra.join(', ') ?: '(none)'}"
                )
            }

            log.info("OK - kept ${kept.size()} record(s): ${kept.join(', ')}")
        }
}
