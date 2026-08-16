nextflow.enable.dsl = 2

include { GUNZIP } from '../../modules/gunzip/main.nf'

process COMPRESS {
    output:
    path "test_file.txt.gz"
    
    script:
    """
    echo "test_content" > test_file.txt
    gzip test_file.txt
    """
}

workflow {
    def meta = [id: 'test_gunzip']
    
    COMPRESS()
    
    GUNZIP(COMPRESS.out.map{ gzip_file -> tuple(meta, gzip_file) })
    
    GUNZIP.out.gunzip.view()
}
