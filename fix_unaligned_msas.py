#!/usr/bin/env python3
"""
Fix unaligned MSA files by re-aligning with MAFFT.

Usage:
    conda activate rna_msa
    python fix_unaligned_msas.py
"""

import subprocess
from pathlib import Path
import tempfile

def is_aligned(a3m_file):
    """Check if all sequences in MSA have same length"""
    sequences = []
    current_seq = ""
    
    with open(a3m_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                current_seq = ""
            else:
                current_seq += line.strip()
        if current_seq:
            sequences.append(current_seq)
    
    if not sequences:
        return True  # Empty file, skip
    
    seq_lengths = set(len(s) for s in sequences)
    return len(seq_lengths) == 1


def realign_with_mafft(a3m_file):
    """Re-align sequences using MAFFT"""
    # Read sequences
    sequences = {}
    current_id = None
    current_seq = ""
    
    with open(a3m_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = current_seq
                current_id = line.strip()[1:].split()[0]
                current_seq = ""
            else:
                current_seq += line.strip()
        if current_id:
            sequences[current_id] = current_seq
    
    if len(sequences) < 2:
        return False, "too_few_sequences"
    
    # Write to temp FASTA (unaligned)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
        for seq_id, seq in sequences.items():
            # Remove any existing gaps
            ungapped = seq.replace('-', '').replace('.', '')
            tmp.write(f">{seq_id}\n{ungapped}\n")
        tmp_path = tmp.name
    
    try:
        # Run MAFFT
        result = subprocess.run(
            ['mafft', '--auto', '--thread', '2', tmp_path],
            capture_output=True,
            text=True,
            timeout=300
        )
        
        if result.returncode != 0:
            return False, "mafft_failed"
        
        # Parse aligned output
        aligned_seqs = {}
        current_id = None
        current_seq = ""
        
        for line in result.stdout.strip().split('\n'):
            if line.startswith('>'):
                if current_id:
                    aligned_seqs[current_id] = current_seq
                current_id = line[1:].split()[0]
                current_seq = ""
            else:
                current_seq += line.strip()
        if current_id:
            aligned_seqs[current_id] = current_seq
        
        # Get query ID (first in original file)
        query_id = list(sequences.keys())[0]
        
        # Write back to A3M with query first
        with open(a3m_file, 'w') as f:
            # Query first
            if query_id in aligned_seqs:
                f.write(f">{query_id}\n{aligned_seqs[query_id].upper()}\n")
            # Then others
            for seq_id, seq in aligned_seqs.items():
                if seq_id != query_id:
                    f.write(f">{seq_id}\n{seq.upper()}\n")
        
        return True, "realigned"
        
    except subprocess.TimeoutExpired:
        return False, "timeout"
    except Exception as e:
        return False, str(e)[:50]
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main():
    msa_dir = Path("RNA-Protein/data/structural_features/msas")
    
    print("Checking for unaligned MSAs...")
    
    unaligned = []
    aligned = 0
    
    for a3m_file in msa_dir.glob("*.a3m"):
        if is_aligned(a3m_file):
            aligned += 1
        else:
            unaligned.append(a3m_file)
    
    print(f"Already aligned: {aligned}")
    print(f"Need realignment: {len(unaligned)}")
    print()
    
    if not unaligned:
        print("Nothing to fix!")
        return
    
    print("Re-aligning with MAFFT...")
    success = 0
    failed = 0
    
    for i, a3m_file in enumerate(unaligned, 1):
        rna_id = a3m_file.stem
        print(f"[{i}/{len(unaligned)}] {rna_id}...", end='', flush=True)
        
        ok, reason = realign_with_mafft(a3m_file)
        
        if ok:
            success += 1
            print(f" ✓")
        else:
            failed += 1
            print(f" ✗ {reason}")
    
    print()
    print(f"Done: {success} fixed, {failed} failed")


if __name__ == "__main__":
    main()
