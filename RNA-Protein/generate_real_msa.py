#!/usr/bin/env python3
"""
Generate a real MSA using BLAST + alignment
"""

import sys
from pathlib import Path
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

output_dir = Path('data/test_real_msa')
output_dir.mkdir(exist_ok=True)

print(f"Generating MSA for: {test_query_id}")
print(f"Output directory: {output_dir}")
print()

# Step 1: BLAST search
print("Step 1: BLAST search for homologs...")
query_file = output_dir / f"{test_query_id}_query.fa"
with open(query_file, 'w') as f:
    f.write(f">{test_query_id}\n{test_query_seq}\n")

blast_out = output_dir / f"{test_query_id}_blast.txt"
cmd = [
    'blastn',
    '-query', str(query_file),
    '-db', 'data/rnacentral/rnacentral_active',
    '-evalue', '0.001',
    '-outfmt', '6 sseqid',
    '-max_target_seqs', '50',  # Get top 50 homologs
    '-num_threads', '4'
]

result = subprocess.run(cmd, capture_output=True, text=True, check=True)
seq_ids = [line.strip() for line in result.stdout.strip().split('\n') if line.strip()]
print(f"  Found {len(seq_ids)} homologs")

# Step 2: Extract sequences
print("Step 2: Extracting homolog sequences...")
id_file = output_dir / f"{test_query_id}_ids.txt"
with open(id_file, 'w') as f:
    for seq_id in seq_ids:
        f.write(f"{seq_id}\n")

homologs_file = output_dir / f"{test_query_id}_homologs.fa"
cmd = [
    'blastdbcmd',
    '-db', 'data/rnacentral/rnacentral_active',
    '-entry_batch', str(id_file),
    '-out', str(homologs_file)
]

subprocess.run(cmd, capture_output=True, text=True, check=True)
print(f"  Extracted to: {homologs_file}")

# Step 3: Align with MAFFT
print("Step 3: Aligning sequences with MAFFT...")
aligned_file = output_dir / f"{test_query_id}_aligned.fa"

cmd = [
    'mafft',
    '--auto',
    '--thread', '4',
    str(homologs_file)
]

result = subprocess.run(cmd, capture_output=True, text=True, check=True)
with open(aligned_file, 'w') as f:
    f.write(result.stdout)

print(f"  Aligned MSA: {aligned_file}")

# Step 4: Convert to A3M format
print("Step 4: Converting to A3M format...")
from src.rna_msa.alignment_formats import fasta_to_a3m, validate_a3m

a3m_file = output_dir / f"{test_query_id}_aligned.a3m"
fasta_to_a3m(aligned_file, a3m_file, query_id=test_query_id)
print(f"  A3M file: {a3m_file}")

# Step 5: Validate and get statistics
print("\nStep 5: Validating A3M...")
stats = validate_a3m(a3m_file)

if stats['valid']:
    print("✓ Valid A3M file")
    print(f"\nMSA Statistics:")
    print(f"  Number of sequences: {stats['n_sequences']}")
    print(f"  Query ID: {stats['query_id']}")
    print(f"  Query length: {stats['query_length']} nt")
    print(f"  Sequence length range: {stats['min_length']}-{stats['max_length']} nt")
    print(f"  Average length: {stats['avg_length']:.1f} nt")
    
    print(f"\n✓ SUCCESS: Generated high-quality MSA!")
    print(f"  MSA depth: {stats['n_sequences']} sequences")
    print(f"  Output: {a3m_file}")
else:
    print(f"✗ Invalid A3M: {stats['error']}")

print("\nFiles generated:")
print(f"  Query: {query_file}")
print(f"  Homologs: {homologs_file}")
print(f"  Aligned FASTA: {aligned_file}")
print(f"  A3M format: {a3m_file}")
