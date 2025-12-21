#!/usr/bin/env python3
"""
Check progress of RNA structure generation

Usage:
    python check_rna_progress.py
"""

import json
from pathlib import Path
from datetime import datetime

def main():
    progress_file = Path("RNA-Protein/data/structural_features/rna_generation_progress.json")
    
    if not progress_file.exists():
        print("No progress file found. Run generate_rna_structures_rpi7317.py first.")
        return
    
    with open(progress_file) as f:
        progress = json.load(f)
    
    completed = progress.get("completed", [])
    failed = progress.get("failed", [])
    last_updated = progress.get("last_updated", "unknown")
    
    total_rnas = 1874
    remaining = total_rnas - len(completed)
    
    print("=" * 60)
    print("RNA Structure Generation Progress")
    print("=" * 60)
    print()
    print(f"Total RNAs: {total_rnas}")
    print(f"Completed: {len(completed)} ({len(completed)/total_rnas*100:.1f}%)")
    print(f"Failed: {len(failed)}")
    print(f"Remaining: {remaining}")
    print()
    print(f"Last updated: {last_updated}")
    print()
    
    # Check actual files
    pdb_dir = Path("RNA-Protein/data/structural_features/rna_structures")
    if pdb_dir.exists():
        pdb_files = list(pdb_dir.glob("*.pdb"))
        print(f"PDB files on disk: {len(pdb_files)}")
    
    msa_dir = Path("RNA-Protein/data/structural_features/msas")
    if msa_dir.exists():
        a3m_files = list(msa_dir.glob("*.a3m"))
        print(f"MSA files on disk: {len(a3m_files)}")
    
    print()
    print("=" * 60)

if __name__ == "__main__":
    main()
