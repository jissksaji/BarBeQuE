nextflow.enable.dsl = 2

include { UNTAR } from '../../modules/untar/main.nf'

process COMPRESS_TAR {
    output:
    path "dummy.tar"
    
    script:
    """
    echo "dummy content" > dummy_test.txt
    tar -cvf dummy.tar dummy_test.txt
    """
}

workflow {
    def meta = [id: 'test_untar']
    
    COMPRESS_TAR()
    
    UNTAR(COMPRESS_TAR.out.map{ tar_file -> tuple(meta, tar_file) })
    
    UNTAR.out.untar.view()
}
