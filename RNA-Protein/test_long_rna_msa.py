#!/usr/bin/env python3
"""
Test MSA generation for a long RNA (>1000nt)
"""

import sys
from pathlib import Path
import subprocess

# Find a long RNA in RNAcentral
print("Finding a long RNA in RNAcentral...")
cmd = "head -1000 data/rnacentral/rnacentral_active.fasta | awk '/^>/ {if (seq) print length(seq), header; header=$0; seq=\"\"} !/^>/ {seq=seq$0} END {print length(seq), header}' | sort -rn | head -5"

result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
print("Longest RNAs found:")
print(result.stdout)

# Let's manually pick a long rRNA (23S rRNA is typically ~2900nt)
# Search for one in the database
print("\nSearching for 23S rRNA (long, conserved RNA)...")
cmd = "grep -m 1 '23S' data/rnacentral/rnacentral_active.fasta"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

if result.returncode == 0:
    header = result.stdout.strip()
    seq_id = header.split()[0][1:]  # Remove '>'
    print(f"Found: {header}")
    print(f"ID: {seq_id}")
    
    # Extract the full sequence
    cmd = f"blastdbcmd -db data/rnacentral/rnacentral_active -entry {seq_id}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    lines = result.stdout.strip().split('\n')
    seq = ''.join(lines[1:])  # Skip header
    
    print(f"Length: {len(seq)} nt")
    
    if len(seq) > 1000:
        print(f"\n✓ This is a long RNA (>{len(seq)} nt)")
        print("\nFor long RNAs, we have two strategies:")
        print("1. Domain splitting: Split into ~250nt chunks and align each")
        print("2. Full-length alignment: If homologs exist, align the whole thing")
        print("\nLet's try full-length alignment first...")
        
        # Test BLAST search
        output_dir = Path('data/test_long_rna_msa')
        output_dir.mkdir(exist_ok=True)
        
        query_file = output_dir / f"{seq_id}_query.fa"
        with open(query_file, 'w') as f:
            f.write(f">{seq_id}\n{seq}\n")
        
        print(f"\nRunning BLAST search...")
        cmd = [
            'blastn',
            '-query', str(query_file),
            '-db', 'data/rnacentral/rnacentral_active',
            '-evalue', '0.001',
            '-outfmt', '6 sseqid pident length',
            '-max_target_seqs', '20',
            '-num_threads', '4'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        hits = [line for line in result.stdout.strip().split('\n') if line]
        
        print(f"Found {len(hits)} homologs")
        print("\nTop 10 hits:")
        print(f"{'Sequence ID':<30} {'Identity':<10} {'Align Length':<15}")
        print("-" * 60)
        for line in hits[:10]:
            parts = line.split('\t')
            print(f"{parts[0]:<30} {parts[1]:>6}%    {parts[2]:>10} nt")
        
        if len(hits) >= 10:
            print(f"\n✓ SUCCESS: Long RNA has enough homologs for MSA!")
            print(f"  Can generate MSA with {len(hits)} sequences")
        else:
            print(f"\n⚠ Limited homologs found")
            print(f"  May need domain splitting strategy")
    else:
        print(f"RNA is only {len(seq)} nt, not long enough for test")
else:
    print("No 23S rRNA found, trying alternative...")
    # Just use a longer sequence from the file
    cmd = "head -10000 data/rnacentral/rnacentral_active.fasta | awk '/^>/ {if (seq && length(seq) > 1000) {print header; print seq; exit}} !/^>/ {if (header) seq=seq$0} /^>/ {header=$0; seq=\"\"}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.stdout:
        lines = result.stdout.strip().split('\n')
        header = lines[0]
        seq = ''.join(lines[1:])
        seq_id = header.split()[0][1:]
        
        print(f"Found long RNA: {header}")
        print(f"Length: {len(seq)} nt")
        print(f"\n✓ Can test long RNA MSA generation with this sequence")
