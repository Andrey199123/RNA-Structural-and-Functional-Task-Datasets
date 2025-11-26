#!/usr/bin/env python3
"""
Test BLAST-based MSA generation on conserved RNAs
"""

import sys
from pathlib import Path
from src.rna_msa.fast_homology_search import FastHomologySearch
import subprocess
import tempfile

# Test query: first rRNA from RNAcentral
test_query_id = "URS000149A9AF"
test_query_seq = """ATACGTAGGTGACAAGCGTTATCCGGATTTACTGGGCGTAAAGGGCGTGTAGGCGGTCTT
GCAAGTCAGAAGTGAAATTCCTGAGCTCAACTCGGGCGCTGCTTCTGAAACTGCAGGACT
TGAGTGCTGGAGGGGATAGCGGAATTCCTAGTGGAGCGGTAAAATGCGCAGATATTAGGA
AGAACACCGGTGGCGAAGGCGGCTATCTGGACAGTAACTGACGCTGAGGCGCGAAATCGT
GGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATACTAGGTG
TAGGGGGTATCGACTCCCTCTGTGCCGCAGTTAACACAATAAGTATTCCGCCTGGGGAGT
ACGGCCGCAAGGTTGAAACTTAAAGGAATTGGCGGGG""".replace('\n', '')

print(f"Testing BLAST search for: {test_query_id}")
print(f"Sequence length: {len(test_query_seq)} nt")
print(f"Type: rRNA (should have many homologs)")
print()

# Create temporary query file
with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
    f.write(f">{test_query_id}\n{test_query_seq}\n")
    query_file = f.name

try:
    # Run BLAST search directly
    print("Running BLAST search...")
    cmd = [
        'blastn',
        '-query', query_file,
        '-db', 'data/rnacentral/rnacentral_active',
        '-evalue', '0.001',
        '-outfmt', '6 sseqid pident length evalue',
        '-max_target_seqs', '100',
        '-num_threads', '4'
    ]
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )
    
    # Parse results
    hits = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) >= 4:
            seq_id = parts[0]
            identity = float(parts[1])
            length = int(parts[2])
            evalue = float(parts[3])
            hits.append((seq_id, identity, length, evalue))
    
    print(f"\n✓ Found {len(hits)} homologs")
    print("\nTop 10 hits:")
    print(f"{'Sequence ID':<30} {'Identity':<10} {'Length':<10} {'E-value':<15}")
    print("-" * 70)
    for seq_id, identity, length, evalue in hits[:10]:
        print(f"{seq_id:<30} {identity:>6.1f}%    {length:>6} nt   {evalue:<15.2e}")
    
    # Statistics
    identities = [h[1] for h in hits]
    print(f"\nIdentity statistics:")
    print(f"  Min: {min(identities):.1f}%")
    print(f"  Max: {max(identities):.1f}%")
    print(f"  Avg: {sum(identities)/len(identities):.1f}%")
    
    high_quality = [h for h in hits if h[1] >= 70]
    print(f"\nHigh-quality homologs (≥70% identity): {len(high_quality)}")
    
    if len(high_quality) >= 10:
        print("\n✓ SUCCESS: Found enough high-quality homologs for MSA!")
        print("  This RNA can generate a good MSA")
    else:
        print("\n⚠ WARNING: Few high-quality homologs")
        print("  MSA quality may be limited")
    
finally:
    Path(query_file).unlink(missing_ok=True)
