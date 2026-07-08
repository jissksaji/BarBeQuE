process test {
    script:
    def format_primer = { seq, f5, f3 ->
        if (seq.contains('#') || seq.length() < (f5 + f3) || (f5 == 0 && f3 == 0)) return seq
        def res = ""
        if (f5 > 0) res += seq[0..<f5].toList().join('#') + '#'
        def mid_start = f5
        def mid_end = seq.length() - f3 - 1
        if (mid_end >= mid_start) res += seq[mid_start..mid_end]
        if (f3 > 0) res += seq[-f3..-1].toList().join('#') + '#'
        return res
    }

    println format_primer("ACTG", 0, 3)
    println format_primer("ACG", 0, 3)
    println format_primer("AC", 0, 3)
    println format_primer("ACTGC", 1, 3)
    println format_primer("ACTG", 0, 0)
    
    """
    echo "done"
    """
}

workflow {
    test()
}
