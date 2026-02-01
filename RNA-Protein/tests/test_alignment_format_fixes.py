#!/usr/bin/env python3
"""
Comprehensive unit tests for alignment format fixes.
Tests every requirement specified by Azam.

Requirements tested:
1. BLAST+MAFFT: Query included in MAFFT input, alignment columns preserved
2. Infernal+Rfam: Stockholm alignment columns preserved
3. Ambiguity codes normalized to 'N'
4. Each sequence on ONE line
5. Query sequence first with NO gaps
6. All sequences have same length
"""

import unittest
import tempfile
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from rna_msa.alignment_formats import (
    fasta_to_a3m, 
    stockholm_to_a3m, 
    validate_a3m,
    parse_stockholm
)


class TestFastaToA3M(unittest.TestCase):
    """Test BLAST+MAFFT branch fixes."""
    
    def setUp(self):
        """Create temporary directory for test files."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
    
    def tearDown(self):
        """Clean up temporary directory."""
        self.temp_dir.cleanup()
    
    def test_query_has_no_gaps(self):
        """Test that query sequence (first) has NO gaps."""
        # Create aligned FASTA with gaps in query
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">query\n")
            f.write("ACGU-ACGU-ACGU\n")  # Query with gaps
            f.write(">seq1\n")
            f.write("ACGUAACGUAACGU\n")
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Read output
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        query_seq = lines[1].strip()
        
        # Query should have NO gaps
        self.assertNotIn('-', query_seq, "Query sequence contains gaps")
        self.assertEqual(query_seq, "ACGUACGUACGU", "Query gaps not removed")
    
    def test_all_sequences_same_length(self):
        """Test that all sequences have the same length (alignment preserved)."""
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">query\n")
            f.write("ACGUACGUACGU\n")
            f.write(">seq1\n")
            f.write("ACGU--GUACGU\n")
            f.write(">seq2\n")
            f.write("----ACGUACGU\n")
            f.write(">seq3\n")
            f.write("ACGUACGU----\n")
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Read sequences
        sequences = {}
        current_id = None
        with open(a3m_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_id = line[1:]
                    sequences[current_id] = ""
                else:
                    sequences[current_id] += line
        
        # All sequences should have same length
        lengths = [len(seq) for seq in sequences.values()]
        self.assertEqual(len(set(lengths)), 1, 
                        f"Sequences have different lengths: {lengths}")
        
        # All should equal query length
        query_len = len(sequences['query'])
        for seq_id, seq in sequences.items():
            self.assertEqual(len(seq), query_len,
                           f"{seq_id} length {len(seq)} != query length {query_len}")
    
    def test_alignment_columns_preserved(self):
        """Test that alignment columns (gaps) are preserved."""
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">query\n")
            f.write("ACGUACGU\n")
            f.write(">seq1\n")
            f.write("AC--ACGU\n")  # 2 gaps in middle
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Read seq1
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        seq1 = lines[3].strip()
        
        # Gaps should be preserved
        self.assertIn('-', seq1, "Gaps were removed from alignment")
        self.assertEqual(seq1.count('-'), 2, "Wrong number of gaps")
        self.assertEqual(len(seq1), 8, "Sequence length changed")
    
    def test_ambiguity_codes_normalized(self):
        """Test that IUPAC ambiguity codes are normalized to 'N'."""
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">query\n")
            f.write("ACGUACGU\n")
            f.write(">seq1\n")
            f.write("WRYKMSBD\n")  # All ambiguity codes
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Read seq1
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        seq1 = lines[3].strip()
        
        # All ambiguity codes should be 'N'
        ambiguity_codes = 'WRYKMSBDHV'
        for code in ambiguity_codes:
            self.assertNotIn(code, seq1, f"Ambiguity code {code} not normalized")
        
        # Should be all N's
        self.assertEqual(seq1, "NNNNNNNN", "Ambiguity codes not normalized to N")
    
    def test_sequences_on_single_lines(self):
        """Test that each sequence is written on ONE line (not wrapped)."""
        # Create long sequence
        long_seq = "ACGU" * 100  # 400 nt
        
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">query\n")
            f.write(long_seq + "\n")
            f.write(">seq1\n")
            f.write(long_seq + "\n")
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Read file
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        # Should have exactly 4 lines: 2 headers + 2 sequences
        self.assertEqual(len(lines), 4, 
                        f"Expected 4 lines, got {len(lines)} (sequences may be wrapped)")
        
        # Check each sequence is on one line
        self.assertTrue(lines[0].startswith('>'), "Line 0 should be header")
        self.assertFalse(lines[1].startswith('>'), "Line 1 should be sequence")
        self.assertTrue(lines[2].startswith('>'), "Line 2 should be header")
        self.assertFalse(lines[3].startswith('>'), "Line 3 should be sequence")
    
    def test_query_is_first(self):
        """Test that query sequence is first in output."""
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">seq1\n")
            f.write("ACGUACGU\n")
            f.write(">query\n")
            f.write("ACGUACGU\n")
            f.write(">seq2\n")
            f.write("ACGUACGU\n")
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Read first header
        with open(a3m_file, 'r') as f:
            first_line = f.readline().strip()
        
        self.assertEqual(first_line, ">query", "Query is not first sequence")
    
    def test_dots_converted_to_dashes(self):
        """Test that Stockholm-style gaps (.) are converted to standard gaps (-)."""
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">query\n")
            f.write("ACGUACGU\n")
            f.write(">seq1\n")
            f.write("AC..ACGU\n")  # Dots as gaps
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Read seq1
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        seq1 = lines[3].strip()
        
        # Dots should be converted to dashes
        self.assertNotIn('.', seq1, "Dots not converted to dashes")
        self.assertIn('-', seq1, "No dashes found")
        self.assertEqual(seq1.count('-'), 2, "Wrong number of gaps")


class TestStockholmToA3M(unittest.TestCase):
    """Test Infernal+Rfam branch fixes."""
    
    def setUp(self):
        """Create temporary directory for test files."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
    
    def tearDown(self):
        """Clean up temporary directory."""
        self.temp_dir.cleanup()
    
    def test_stockholm_alignment_preserved(self):
        """Test that Stockholm alignment columns are preserved."""
        # Create Stockholm file with alignment
        sto_file = self.temp_path / "test.sto"
        with open(sto_file, 'w') as f:
            f.write("# STOCKHOLM 1.0\n")
            f.write("#=GF AC RF00001\n")
            f.write("seq1    ACGU-ACGU\n")
            f.write("seq2    AC--AACGU\n")
            f.write("seq3    ACGUAAC--\n")
            f.write("//\n")
        
        a3m_file = self.temp_path / "output.a3m"
        stockholm_to_a3m(sto_file, a3m_file, query_id="seq1")
        
        # Read sequences
        sequences = {}
        current_id = None
        with open(a3m_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_id = line[1:]
                    sequences[current_id] = ""
                else:
                    sequences[current_id] += line
        
        # All sequences should have same length
        lengths = [len(seq) for seq in sequences.values()]
        self.assertEqual(len(set(lengths)), 1,
                        f"Stockholm alignment not preserved: {lengths}")
    
    def test_stockholm_query_no_gaps(self):
        """Test that query from Stockholm has no gaps."""
        sto_file = self.temp_path / "test.sto"
        with open(sto_file, 'w') as f:
            f.write("# STOCKHOLM 1.0\n")
            f.write("query    ACGU-ACGU\n")
            f.write("seq1     ACGUAACGU\n")
            f.write("//\n")
        
        a3m_file = self.temp_path / "output.a3m"
        stockholm_to_a3m(sto_file, a3m_file, query_id="query")
        
        # Read query
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        query_seq = lines[1].strip()
        
        # Query should have no gaps
        self.assertNotIn('-', query_seq, "Query has gaps")
        self.assertEqual(query_seq, "ACGUACGU", "Query gaps not removed")
    
    def test_stockholm_dots_to_dashes(self):
        """Test that Stockholm gaps (.) are converted to standard gaps (-)."""
        sto_file = self.temp_path / "test.sto"
        with open(sto_file, 'w') as f:
            f.write("# STOCKHOLM 1.0\n")
            f.write("query    ACGUACGU\n")
            f.write("seq1     AC..ACGU\n")  # Dots as gaps
            f.write("//\n")
        
        a3m_file = self.temp_path / "output.a3m"
        stockholm_to_a3m(sto_file, a3m_file, query_id="query")
        
        # Read seq1
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        seq1 = lines[3].strip()
        
        # Dots should be dashes
        self.assertNotIn('.', seq1, "Dots not converted")
        self.assertIn('-', seq1, "No dashes found")
    
    def test_stockholm_ambiguity_normalized(self):
        """Test that ambiguity codes in Stockholm are normalized."""
        sto_file = self.temp_path / "test.sto"
        with open(sto_file, 'w') as f:
            f.write("# STOCKHOLM 1.0\n")
            f.write("query    ACGUACGU\n")
            f.write("seq1     WRYKMSBD\n")  # Ambiguity codes
            f.write("//\n")
        
        a3m_file = self.temp_path / "output.a3m"
        stockholm_to_a3m(sto_file, a3m_file, query_id="query")
        
        # Read seq1
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        seq1 = lines[3].strip()
        
        # Should be all N's
        ambiguity_codes = 'WRYKMSBDHV'
        for code in ambiguity_codes:
            self.assertNotIn(code, seq1, f"Ambiguity code {code} not normalized")
    
    def test_stockholm_external_query(self):
        """Test Stockholm with external query sequence (Rfam case)."""
        sto_file = self.temp_path / "test.sto"
        with open(sto_file, 'w') as f:
            f.write("# STOCKHOLM 1.0\n")
            f.write("seed1    ACGUACGU\n")
            f.write("seed2    AC--ACGU\n")
            f.write("//\n")
        
        # External query not in seed alignment
        query_seq = "ACGUACGU"
        
        a3m_file = self.temp_path / "output.a3m"
        stockholm_to_a3m(sto_file, a3m_file, 
                        query_id="my_query", 
                        query_sequence=query_seq)
        
        # Read first sequence
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        # First should be query
        self.assertEqual(lines[0].strip(), ">my_query", "Query not first")
        query_out = lines[1].strip()
        self.assertEqual(query_out, query_seq, "Query sequence wrong")
        self.assertNotIn('-', query_out, "Query has gaps")


class TestValidateA3M(unittest.TestCase):
    """Test validation function."""
    
    def setUp(self):
        """Create temporary directory for test files."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
    
    def tearDown(self):
        """Clean up temporary directory."""
        self.temp_dir.cleanup()
    
    def test_valid_a3m(self):
        """Test that valid A3M passes validation."""
        a3m_file = self.temp_path / "valid.a3m"
        with open(a3m_file, 'w') as f:
            f.write(">query\n")
            f.write("ACGUACGU\n")
            f.write(">seq1\n")
            f.write("AC--ACGU\n")
            f.write(">seq2\n")
            f.write("ACGUAC--\n")
        
        stats = validate_a3m(a3m_file)
        
        self.assertTrue(stats['valid'], f"Valid A3M failed: {stats.get('error')}")
        self.assertEqual(stats['n_sequences'], 3)
        self.assertEqual(stats['query_length'], 8)
    
    def test_invalid_query_with_gaps(self):
        """Test that query with gaps fails validation."""
        a3m_file = self.temp_path / "invalid.a3m"
        with open(a3m_file, 'w') as f:
            f.write(">query\n")
            f.write("AC--ACGU\n")  # Query has gaps
            f.write(">seq1\n")
            f.write("ACGUACGU\n")
        
        stats = validate_a3m(a3m_file)
        
        self.assertFalse(stats['valid'], "Query with gaps should fail")
        self.assertIn('gaps', stats['error'].lower())
    
    def test_invalid_different_lengths(self):
        """Test that different sequence lengths fail validation."""
        a3m_file = self.temp_path / "invalid.a3m"
        with open(a3m_file, 'w') as f:
            f.write(">query\n")
            f.write("ACGUACGU\n")  # Length 8
            f.write(">seq1\n")
            f.write("ACGU\n")  # Length 4 - WRONG
        
        stats = validate_a3m(a3m_file)
        
        self.assertFalse(stats['valid'], "Different lengths should fail")
        self.assertIn('length', stats['error'].lower())
    
    def test_invalid_ambiguity_codes(self):
        """Test that ambiguity codes fail validation."""
        a3m_file = self.temp_path / "invalid.a3m"
        with open(a3m_file, 'w') as f:
            f.write(">query\n")
            f.write("ACGUACGU\n")
            f.write(">seq1\n")
            f.write("WRYKMSBD\n")  # Ambiguity codes
        
        stats = validate_a3m(a3m_file)
        
        self.assertFalse(stats['valid'], "Ambiguity codes should fail")
        self.assertIn('ambiguity', stats['error'].lower())


class TestTrRosettaRNA2Compatibility(unittest.TestCase):
    """Test compatibility with trRosettaRNA2 example format."""
    
    def test_example_format(self):
        """Test that trRosettaRNA2 example passes validation."""
        example_file = Path('trRosettaRNA2/example/seq.a3m')
        
        if not example_file.exists():
            self.skipTest("trRosettaRNA2 example not found")
        
        stats = validate_a3m(example_file)
        
        self.assertTrue(stats['valid'], 
                       f"trRosettaRNA2 example failed: {stats.get('error')}")
    
    def test_generated_matches_example_format(self):
        """Test that our generated A3M matches trRosettaRNA2 format."""
        temp_dir = tempfile.TemporaryDirectory()
        temp_path = Path(temp_dir.name)
        
        try:
            # Create aligned FASTA
            aligned_fasta = temp_path / "aligned.fa"
            with open(aligned_fasta, 'w') as f:
                f.write(">query\n")
                f.write("CCGCCGCGCCAUGCCUGUGGCGGAAACCGCCGCGCCAUGCCUGUGGCGG\n")
                f.write(">seq1\n")
                f.write("------------GCCCGTGGCGGAAACCGCCGCGCC-------------\n")
                f.write(">seq2\n")
                f.write("CCGCCGCGCCGTGCCTGTG-----------CGCGCCGTGCCTGTGGCGC\n")
            
            # Convert to A3M
            a3m_file = temp_path / "output.a3m"
            fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
            
            # Validate
            stats = validate_a3m(a3m_file)
            
            self.assertTrue(stats['valid'], f"Generated A3M invalid: {stats.get('error')}")
            
            # Check specific properties
            self.assertEqual(stats['query_length'], 49, "Wrong query length")
            self.assertTrue(stats['all_same_length'], "Sequences have different lengths")
            
        finally:
            temp_dir.cleanup()


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""
    
    def setUp(self):
        """Create temporary directory for test files."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
    
    def tearDown(self):
        """Clean up temporary directory."""
        self.temp_dir.cleanup()
    
    def test_empty_file(self):
        """Test handling of empty file."""
        empty_file = self.temp_path / "empty.fa"
        empty_file.touch()
        
        a3m_file = self.temp_path / "output.a3m"
        
        with self.assertRaises(ValueError):
            fasta_to_a3m(empty_file, a3m_file)
    
    def test_query_not_found(self):
        """Test error when query ID not in alignment."""
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">seq1\n")
            f.write("ACGUACGU\n")
        
        a3m_file = self.temp_path / "output.a3m"
        
        with self.assertRaises(ValueError):
            fasta_to_a3m(aligned_fasta, a3m_file, query_id="nonexistent")
    
    def test_very_long_sequence(self):
        """Test handling of very long sequences (no wrapping)."""
        long_seq = "ACGU" * 1000  # 4000 nt
        
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">query\n")
            f.write(long_seq + "\n")
            f.write(">seq1\n")
            f.write(long_seq + "\n")
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Check no wrapping
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        self.assertEqual(len(lines), 4, "Sequences were wrapped")
        self.assertEqual(len(lines[1].strip()), 4000, "Query length wrong")
    
    def test_mixed_case(self):
        """Test that mixed case is handled correctly."""
        aligned_fasta = self.temp_path / "aligned.fa"
        with open(aligned_fasta, 'w') as f:
            f.write(">query\n")
            f.write("AcGuAcGu\n")  # Mixed case
            f.write(">seq1\n")
            f.write("acguacgu\n")  # Lowercase
        
        a3m_file = self.temp_path / "output.a3m"
        fasta_to_a3m(aligned_fasta, a3m_file, query_id="query")
        
        # Read sequences
        with open(a3m_file, 'r') as f:
            lines = f.readlines()
        
        query_seq = lines[1].strip()
        seq1 = lines[3].strip()
        
        # Should be uppercase
        self.assertEqual(query_seq, "ACGUACGU", "Not uppercase")
        self.assertEqual(seq1, "ACGUACGU", "Not uppercase")


def run_tests():
    """Run all tests and print summary."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestFastaToA3M))
    suite.addTests(loader.loadTestsFromTestCase(TestStockholmToA3M))
    suite.addTests(loader.loadTestsFromTestCase(TestValidateA3M))
    suite.addTests(loader.loadTestsFromTestCase(TestTrRosettaRNA2Compatibility))
    suite.addTests(loader.loadTestsFromTestCase(TestEdgeCases))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Successes: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    
    if result.wasSuccessful():
        print("\n✅ ALL TESTS PASSED - Ready for production!")
    else:
        print("\n❌ SOME TESTS FAILED - Review failures above")
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
