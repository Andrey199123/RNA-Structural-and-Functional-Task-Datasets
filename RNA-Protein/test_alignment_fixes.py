#!/usr/bin/env python3
"""
Test script to verify alignment format fixes.
Compares output format with trRosettaRNA2 example.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))
from rna_msa.alignment_formats import validate_a3m


def test_a3m_format(a3m_file: Path):
    """Test A3M file format against trRosettaRNA2 requirements."""
    
    print(f"\nTesting: {a3m_file}")
    print("=" * 60)
    
    # Read sequences
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(a3m_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    if not sequences:
        print("❌ FAIL: No sequences found")
        return False
    
    # Get query (first sequence)
    query_id = list(sequences.keys())[0]
    query_seq = sequences[query_id]
    
    print(f"Query ID: {query_id}")
    print(f"Query length: {len(query_seq)}")
    print(f"Total sequences: {len(sequences)}")
    
    # Check 1: Query has no gaps
    if '-' in query_seq:
        print(f"❌ FAIL: Query contains {query_seq.count('-')} gaps")
        return False
    else:
        print("✓ Query has no gaps")
    
    # Check 2: All sequences same length
    seq_lengths = [len(seq) for seq in sequences.values()]
    if len(set(seq_lengths)) > 1:
        print(f"❌ FAIL: Sequences have different lengths: {set(seq_lengths)}")
        print(f"   Min: {min(seq_lengths)}, Max: {max(seq_lengths)}")
        return False
    else:
        print(f"✓ All sequences have same length ({seq_lengths[0]})")
    
    # Check 3: No ambiguity codes
    ambiguity_codes = 'WRYKMSBDHV'
    sequences_with_ambiguity = []
    for seq_id, seq in sequences.items():
        if any(char in ambiguity_codes for char in seq.upper()):
            sequences_with_ambiguity.append(seq_id)
    
    if sequences_with_ambiguity:
        print(f"❌ FAIL: Found ambiguity codes in {len(sequences_with_ambiguity)} sequences")
        print(f"   Examples: {sequences_with_ambiguity[:3]}")
        return False
    else:
        print("✓ No ambiguity codes found")
    
    # Check 4: Sequences on single lines (check file directly)
    with open(a3m_file, 'r') as f:
        lines = f.readlines()
    
    multi_line_seqs = 0
    in_seq = False
    for line in lines:
        if line.startswith('>'):
            in_seq = True
        elif in_seq and line.strip():
            # Check if next line is also sequence (not header)
            idx = lines.index(line)
            if idx + 1 < len(lines):
                next_line = lines[idx + 1].strip()
                if next_line and not next_line.startswith('>'):
                    multi_line_seqs += 1
            in_seq = False
    
    if multi_line_seqs > 0:
        print(f"⚠️  WARNING: Found {multi_line_seqs} sequences split across multiple lines")
    else:
        print("✓ All sequences on single lines")
    
    # Check 5: Gap character consistency
    gap_chars = set()
    for seq in sequences.values():
        for char in seq:
            if char in ['-', '.', '_']:
                gap_chars.add(char)
    
    if len(gap_chars) > 1:
        print(f"⚠️  WARNING: Multiple gap characters found: {gap_chars}")
    elif gap_chars == {'-'}:
        print("✓ Consistent gap character (-)")
    elif not gap_chars:
        print("✓ No gaps in alignment")
    
    # Check 6: Alignment structure
    gap_counts = [seq.count('-') for seq in sequences.values()]
    if len(sequences) > 1:
        avg_gaps = sum(gap_counts[1:]) / (len(gap_counts) - 1)  # Exclude query
        print(f"✓ Average gaps per sequence: {avg_gaps:.1f}")
    
    # Validate using built-in function
    print("\nRunning validate_a3m()...")
    stats = validate_a3m(a3m_file)
    
    if stats['valid']:
        print("✅ PASS: A3M file is valid for trRosettaRNA2")
        return True
    else:
        print(f"❌ FAIL: {stats.get('error', 'Unknown error')}")
        return False


def compare_with_example():
    """Compare format with trRosettaRNA2 example."""
    
    example_file = Path('trRosettaRNA2/example/seq.a3m')
    
    if not example_file.exists():
        print(f"Example file not found: {example_file}")
        return
    
    print("\n" + "=" * 60)
    print("COMPARING WITH trRosettaRNA2 EXAMPLE")
    print("=" * 60)
    
    test_a3m_format(example_file)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Test A3M format compatibility")
    parser.add_argument('a3m_files', nargs='*', type=Path, help='A3M files to test')
    parser.add_argument('--example', action='store_true', help='Test trRosettaRNA2 example')
    
    args = parser.parse_args()
    
    if args.example:
        compare_with_example()
    
    if args.a3m_files:
        all_pass = True
        for a3m_file in args.a3m_files:
            if not a3m_file.exists():
                print(f"\n❌ File not found: {a3m_file}")
                all_pass = False
                continue
            
            passed = test_a3m_format(a3m_file)
            if not passed:
                all_pass = False
        
        print("\n" + "=" * 60)
        if all_pass:
            print("✅ ALL TESTS PASSED")
        else:
            print("❌ SOME TESTS FAILED")
        print("=" * 60)
    
    if not args.example and not args.a3m_files:
        print("Usage:")
        print("  python test_alignment_fixes.py file1.a3m file2.a3m")
        print("  python test_alignment_fixes.py --example")
