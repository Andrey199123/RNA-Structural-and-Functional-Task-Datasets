#!/usr/bin/env python3
"""
Download pre-computed protein structures from AlphaFold Database

For RPI7317 proteins with UniProt IDs, download their AlphaFold structures
instead of predicting them. This is much faster and often more accurate.

Usage:
    python download_alphafold_structures.py --mode test  # 10 proteins
    python download_alphafold_structures.py --mode full  # All 118 proteins
"""

import argparse
import time
from pathlib import Path
import requests
from Bio import SeqIO


def download_alphafold_structure(uniprot_id: str, output_dir: Path) -> Path:
    """Download AlphaFold structure for a UniProt ID"""
    output_file = output_dir / f"{uniprot_id}.pdb"
    
    if output_file.exists():
        print(f"  Structure already exists: {uniprot_id}")
        return output_file
    
    print(f"  Downloading AlphaFold structure: {uniprot_id}")
    
    # Try different AlphaFold DB versions (v6 is latest as of 2025)
    versions = ["v6", "v4", "v3", "v2"]
    
    for version in versions:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{version}.pdb"
        
        try:
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                output_file.write_text(response.text)
                print(f"    ✓ Downloaded from AlphaFold DB ({version})")
                return output_file
                
        except Exception:
            continue
    
    print(f"    ✗ Not found in AlphaFold DB")
    return None


def main():
    parser = argparse.ArgumentParser(description="Download AlphaFold structures for RPI7317")
    parser.add_argument("--mode", choices=["test", "full"], default="test",
                       help="test: process 10 proteins, full: process all 118")
    parser.add_argument("--output-dir", type=Path, 
                       default=Path("RNA-Protein/data/structural_features/protein_structures"),
                       help="Output directory for PDB files")
    args = parser.parse_args()
    
    print("=" * 80)
    print("AlphaFold Structure Download for RPI7317")
    print("=" * 80)
    print()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load protein sequences to get UniProt IDs
    print("Loading protein sequences...")
    protein_fasta = Path("RNA-Protein/data/RPI7317/RPI7317_protein_sequence.fasta")
    
    proteins = {}
    for record in SeqIO.parse(protein_fasta, "fasta"):
        # Extract UniProt ID from header (format: >sp|Q01844|EWS_HUMAN ...)
        if "|" in record.id:
            uniprot_id = record.id.split("|")[1]
        else:
            uniprot_id = record.id
        proteins[uniprot_id] = str(record.seq)
    
    print(f"  Found {len(proteins)} proteins")
    print()
    
    # Limit to test samples if in test mode
    if args.mode == "test":
        proteins = dict(list(proteins.items())[:10])
        print(f"TEST MODE: Processing {len(proteins)} proteins")
        print()
    
    # Download structures
    stats = {"success": 0, "not_found": 0, "failed": 0}
    
    for i, (uniprot_id, sequence) in enumerate(proteins.items(), 1):
        print(f"[{i}/{len(proteins)}] Processing: {uniprot_id} ({len(sequence)} aa)")
        
        pdb_file = download_alphafold_structure(uniprot_id, args.output_dir)
        
        if pdb_file:
            stats["success"] += 1
        else:
            # Check if it was not found vs other error
            if "Not found" in str(pdb_file):
                stats["not_found"] += 1
            else:
                stats["failed"] += 1
        
        time.sleep(0.5)  # Rate limiting
        print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print(f"Proteins processed: {len(proteins)}")
    print(f"  - Downloaded from AlphaFold DB: {stats['success']}")
    print(f"  - Not found in AlphaFold DB: {stats['not_found']}")
    print(f"  - Download failed: {stats['failed']}")
    print()
    print(f"Success rate: {stats['success']/len(proteins)*100:.1f}%")
    print()
    print(f"Output directory: {args.output_dir}")
    print()
    
    if stats['not_found'] > 0:
        print("NOTE: Proteins not found in AlphaFold DB will need alternative prediction")
        print("      (e.g., ESMFold local installation or trRosetta2)")
    print()


if __name__ == "__main__":
    main()
