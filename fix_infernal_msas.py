#!/usr/bin/env python3
"""
Fix MSA files that are missing the query sequence.
These are Infernal-generated MSAs where the query wasn't added.

Usage:
    conda activate rna_msa
    python fix_infernal_msas.py
"""

from pathlib import Path
from Bio import SeqIO

def main():
    msa_dir = Path("RNA-Protein/data/structural_features/msas")
    fasta_file = Path("RNA-Protein/data/RPI7317/RPI7317_IncRNA_sequence.fasta")
    
    # Load all RNA sequences
    print("Loading RNA sequences...")
    rna_seqs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        rna_seqs[record.id] = str(record.seq)
    print(f"Loaded {len(rna_seqs)} sequences")
    
    # Find and fix broken MSAs
    fixed = 0
    already_ok = 0
    not_found = 0
    
    for a3m_file in msa_dir.glob("*.a3m"):
        rna_id = a3m_file.stem
        
        # Read first line
        with open(a3m_file) as f:
            first_line = f.readline().strip()
        
        # Check if query is present
        if first_line == f">{rna_id}":
            already_ok += 1
            continue
        
        # Need to fix - prepend query sequence
        if rna_id not in rna_seqs:
            print(f"  Warning: {rna_id} not found in FASTA")
            not_found += 1
            continue
        
        query_seq = rna_seqs[rna_id]
        
        # Read existing content
        with open(a3m_file) as f:
            existing_content = f.read()
        
        # Write with query prepended
        with open(a3m_file, 'w') as f:
            f.write(f">{rna_id}\n{query_seq}\n")
            f.write(existing_content)
        
        print(f"  Fixed: {rna_id}")
        fixed += 1
    
    print()
    print(f"Summary:")
    print(f"  Already OK: {already_ok}")
    print(f"  Fixed: {fixed}")
    print(f"  Not found in FASTA: {not_found}")

if __name__ == "__main__":
    main()
