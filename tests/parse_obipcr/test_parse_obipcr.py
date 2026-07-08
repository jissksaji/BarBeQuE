import unittest
import sys
import os

# Add the bin directory to path so we can import parse_obipcr.py
bin_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../bin'))
sys.path.insert(0, bin_dir)

import parse_obipcr

class TestParseObiPCR(unittest.TestCase):

    def test_parse_fasta_header(self):
        header = '>seq1_sub[10..20] {"forward_primer":"ACGT", "direction":"F", "forward_match":"ACGT"}'
        parsed = parse_obipcr.parse_fasta_header(header)
        self.assertIsNotNone(parsed)
        self.assertEqual(parsed['seq_id'], 'seq1')
        self.assertEqual(parsed['start_pos'], 10)
        self.assertEqual(parsed['end_pos'], 20)
        self.assertEqual(parsed['metadata']['forward_primer'], 'ACGT')
        
        # Invalid header
        self.assertIsNone(parse_obipcr.parse_fasta_header('invalid_header'))

    def test_calculate_gc(self):
        self.assertEqual(parse_obipcr.calculate_gc(''), "0.0")
        self.assertEqual(parse_obipcr.calculate_gc('GCGC'), "100.0")
        self.assertEqual(parse_obipcr.calculate_gc('ATAT'), "0.0")
        self.assertEqual(parse_obipcr.calculate_gc('GCTA'), "50.0")

    def test_is_iupac_match(self):
        # Exact matches
        self.assertTrue(parse_obipcr.is_iupac_match('A', 'A'))
        self.assertTrue(parse_obipcr.is_iupac_match('C', 'C'))
        
        # Degenerate matches
        self.assertTrue(parse_obipcr.is_iupac_match('Y', 'C'))
        self.assertTrue(parse_obipcr.is_iupac_match('Y', 'T'))
        self.assertTrue(parse_obipcr.is_iupac_match('R', 'A'))
        self.assertTrue(parse_obipcr.is_iupac_match('R', 'G'))
        self.assertTrue(parse_obipcr.is_iupac_match('N', 'A'))
        
        # Mismatches
        self.assertFalse(parse_obipcr.is_iupac_match('Y', 'A'))
        self.assertFalse(parse_obipcr.is_iupac_match('R', 'C'))
        self.assertFalse(parse_obipcr.is_iupac_match('A', 'C'))

    def test_calculate_tm(self):
        # Empty sequence
        self.assertEqual(parse_obipcr.calculate_tm(''), "0.0")
        # Short (<14 nt) uses Wallace rule: 2*(A+T) + 4*(G+C)
        self.assertEqual(parse_obipcr.calculate_tm('GCGC'), "16.0")
        self.assertEqual(parse_obipcr.calculate_tm('ACGT'), "12.0")
        # Long (>=14 nt) uses GC formula: 64.9 + 41*(GC-16.4)/N
        self.assertEqual(parse_obipcr.calculate_tm('GCGCGCGCGCGCGCGCGCGC'), "72.3")

    def test_gc_clamp_count(self):
        self.assertEqual(parse_obipcr.gc_clamp_count('A'), "0")
        self.assertEqual(parse_obipcr.gc_clamp_count('AAAAAGCGCC'), "5") # last 5 is GCGCC -> 5 G/C
        self.assertEqual(parse_obipcr.gc_clamp_count('AAAAAATATA'), "0") # last 5 is ATATA -> 0 G/C

    def test_max_mononucleotide_run(self):
        self.assertEqual(parse_obipcr.max_mononucleotide_run(''), "0")
        self.assertEqual(parse_obipcr.max_mononucleotide_run('AGCT'), "1")
        self.assertEqual(parse_obipcr.max_mononucleotide_run('AAGGGCTTTT'), "4")

    def test_reverse_complement(self):
        self.assertEqual(parse_obipcr.reverse_complement(''), "")
        self.assertEqual(parse_obipcr.reverse_complement('ACGTN'), "NACGT")
        self.assertEqual(parse_obipcr.reverse_complement('agct'), "AGCT")

    def test_three_prime_dimer_length(self):
        self.assertEqual(parse_obipcr.three_prime_dimer_length('ACGT', 'ACGT'), "4")
        self.assertEqual(parse_obipcr.three_prime_dimer_length('AAAA', 'TTTT'), "4")
        self.assertEqual(parse_obipcr.three_prime_dimer_length('ACGT', 'AAAA'), "1")  # T matches A
        self.assertEqual(parse_obipcr.three_prime_dimer_length('', 'ACGT'), "0")
        
    def test_max_hairpin_stem(self):
        self.assertEqual(parse_obipcr.max_hairpin_stem('ACGT'), "0")
        # AGCT ... AGCT -> forms stem of 4 with a loop of 3 (TTT)
        self.assertEqual(parse_obipcr.max_hairpin_stem('AGCTTTTAGCT'), "4")
        # No loop (too short)
        self.assertEqual(parse_obipcr.max_hairpin_stem('AAAAA'), "0")
    def test_find_mismatches(self):
        pm, gm, tpm, sev = parse_obipcr.find_mismatches('ACGT', 'ACTT', 100, is_reverse=False)
        self.assertEqual(pm, "3")
        self.assertEqual(gm, "102")
        self.assertEqual(tpm, "2")
        self.assertEqual(sev, "3.0")
        
        # Test IUPAC degenerate bases (Y matches C or T)
        # ACYT vs ACCT (no mismatch)
        pm, gm, tpm, sev = parse_obipcr.find_mismatches('ACYT', 'ACCT', 100, is_reverse=False)
        self.assertEqual(pm, "")
        self.assertEqual(gm, "")
        self.assertEqual(tpm, "")
        self.assertEqual(sev, "0.0")
        
        # Reverse case: is_reverse=True
        pm, gm, tpm, sev = parse_obipcr.find_mismatches('ACGT', 'ACTT', 100, is_reverse=True)
        self.assertEqual(pm, "3")
        self.assertEqual(gm, "98") # 100 - 2
        self.assertEqual(tpm, "2")
        self.assertEqual(sev, "3.0")

    def test_process_obipcr(self):
        import tempfile
        
        fasta_content = '>seq1_sub[100..200] {"forward_primer":"ACGT", "forward_match":"ACGT", "reverse_primer":"TGCA", "reverse_match":"TGCA", "direction":"F", "forward_error": 0, "reverse_error": 0}\n'
        fasta_content += 'A' * 101 + '\n'
        
        with tempfile.NamedTemporaryFile('w', delete=False) as f_in:
            f_in.write(fasta_content)
            in_name = f_in.name
            
        with tempfile.NamedTemporaryFile('w', delete=False) as f_out:
            out_name = f_out.name
            
        try:
            parse_obipcr.process_obipcr(in_name, out_name)
            with open(out_name, 'r') as f:
                lines = f.readlines()
                self.assertEqual(len(lines), 2)
                self.assertTrue(lines[0].startswith('Sequence_ID'))
                self.assertTrue(lines[1].startswith('seq1'))
        finally:
            os.remove(in_name)
            os.remove(out_name)

if __name__ == '__main__':
    unittest.main()
