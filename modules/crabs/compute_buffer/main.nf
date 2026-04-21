//process COMPUTE_BUFFER {
//
//    input:
//    path(db)
//
//    output:
//    env(buffersize)
//
//    script:
//    """
//    buffersize=\$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq\$0} END{print length(seq)*2}' $db | sort -n | tail -1)
//    """
//}