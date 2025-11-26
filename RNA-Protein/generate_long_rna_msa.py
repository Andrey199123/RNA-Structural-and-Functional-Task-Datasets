#!/usr/bin/env python3
"""
Generate MSA for a long RNA (1456 nt)
"""

import subprocess
from pathlib import Path

seq_id = "URS00014B6922"
print(f"Generating MSA for long RNA: {seq_id}")

# Extract sequence
cmd = f"blastdbcmd -db data/rnacentral/rnacentral_active -entry {seq_id}"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
lines = result.stdout.strip().split('\n')
seq = ''.join(lines[1:])

print(f"Length: {len(seq)} nt")
print()

output_dir = Path('data/test_long_rna_msa')
output_dir.mkdir(exist_ok=True)

# Step 1: BLAST search
print("Step 1: BLAST search...")
query_file = output_dir / f"{seq_id}_query.fa"
with open(query_file, 'w') as f:
    f.write(f">{seq_id}\n{seq}\n")

cmd = [
    'blastn',
    '-query', str(query_file),
    '-db', 'data/rnacentral/rnacentral_active',
    '-evalue', '0.001',
    '-outfmt', '6 sseqid pident',
    '-max_target_seqs', '30',
    '-num_threads', '4'
]

result = subprocess.run(cmd, capture_output=True, text=True, check=True)
hits = []
for line in result.stdout.strip().split('\n'):
    if line:
        parts = line.split('\t')
        seq_id_hit = parts[0]
        identity = float(parts[1])
        if identity >= 70:  # Only high-quality homologs
            hits.append(seq_id_hit)

print(f"  Found {len(hits)} high-quality homologs (≥70% identity)")

if len(hits) < 5:
    print("\n⚠ Too few homologs for good MSA")
    print("  Would use domain splitting strategy instead")
    exit(0)

# Step 2: Extract sequences
print("Step 2: Extracting sequences...")
id_file = output_dir / f"{seq_id}_ids.txt"
with open(id_file, 'w') as f:
    for hit_id in hits:
        f.write(f"{hit_id}\n")

homologs_file = output_dir / f"{seq_id}_homologs.fa"
cmd = [
    'blastdbcmd',
    '-db', 'data/rnacentral/rnacentral_active',
    '-entry_batch', str(id_file),
    '-out', str(homologs_file)
]
subprocess.run(cmd, capture_output=True, text=True, check=True)

# Step 3: Align with MAFFT
print("Step 3: Aligning with MAFFT (this may take a minute for long sequences)...")
aligned_file = output_dir / f"{seq_id}_aligned.fa"

cmd = [
    'mafft',
    '--auto',
    '--thread', '4',
    str(homologs_file)
]

result = subprocess.run(cmd, capture_output=True, text=True, check=True)
with open(aligned_file, 'w') as f:
    f.write(result.stdout)

# Step 4: Convert to A3M
print("Step 4: Converting to A3M...")
from src.rna_msa.alignment_formats import fasta_to_a3m, validate_a3m

a3m_file = output_dir / f"{seq_id}_aligned.a3m"
fasta_to_a3m(aligned_file, a3m_file, query_id=seq_id)

# Step 5: Validate
print("\nStep 5: Validation...")
stats = validate_a3m(a3m_file)

if stats['valid']:
    print("✓ Valid A3M file")
    print(f"\nMSA Statistics:")
    print(f"  Number of sequences: {stats['n_sequences']}")
    print(f"  Query length: {stats['query_length']} nt")
    print(f"  Length range: {stats['min_length']}-{stats['max_length']} nt")
    
    print(f"\n✓ SUCCESS: Generated MSA for long RNA!")
    print(f"  Length: {stats['query_length']} nt (>1000nt)")
    print(f"  MSA depth: {stats['n_sequences']} sequences")
    print(f"  Output: {a3m_file}")
else:
    print(f"✗ Invalid: {stats['error']}")
