nextflow.enable.dsl = 2

include { CAT_FASTQ } from '../../modules/cat_fastq/main.nf'

process GENERATE_GZ {
    output:
    path "dummy_*.fastq.gz"

    script:
    """
    echo "@SEQ" > dummy_1.fastq
    echo "GATC" >> dummy_1.fastq
    echo "+" >> dummy_1.fastq
    echo "!!!!" >> dummy_1.fastq
    gzip -c dummy_1.fastq > dummy_1.fastq.gz
    
    echo "@SEQ" > dummy_2.fastq
    echo "GATC" >> dummy_2.fastq
    echo "+" >> dummy_2.fastq
    echo "!!!!" >> dummy_2.fastq
    gzip -c dummy_2.fastq > dummy_2.fastq.gz
    """
}

workflow {
    def meta = [sample_id: 'test_sample']
    
    GENERATE_GZ()
    
    // We need a list of reads
    def reads = GENERATE_GZ.out.flatten().collate(2)
    
    CAT_FASTQ(reads.map{ r -> tuple(meta, r) })
    
    CAT_FASTQ.out.reads.view()
}
