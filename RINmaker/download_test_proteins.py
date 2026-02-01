#!/usr/bin/env python3
"""
Download test protein structures mentioned in RINmaker paper for benchmarking
"""

import os
import urllib.request
from pathlib import Path

# Test structures from RINmaker paper (Table 2)
TEST_PROTEINS = {
    '1al1': {'residues': 12, 'description': 'Very small protein'},
    '1l2y': {'residues': 20, 'description': 'Trp-cage miniprotein'},
    '1bk8': {'residues': 50, 'description': 'Small protein'},
    '1eod': {'residues': 100, 'description': 'Medium protein'},
    '1hhq': {'residues': 150, 'description': 'Medium protein'},
    '1ab5': {'residues': 250, 'description': 'Large protein'},
    '1bge': {'residues': 318, 'description': 'Large protein'},
    '1b1y': {'residues': 500, 'description': 'Very large protein'},
    '2vji': {'residues': 598, 'description': 'Very large protein'},
    '1miq': {'residues': 750, 'description': 'Very large protein'},
    '1a3w': {'residues': 1000, 'description': 'Extra large protein'},
}

def download_pdb(pdb_id, info, output_dir='test_proteins'):
    """Download PDB file from RCSB"""
    Path(output_dir).mkdir(exist_ok=True)
    
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    output_file = os.path.join(output_dir, f'{pdb_id}.pdb')
    
    if os.path.exists(output_file):
        print(f"✓ {pdb_id} already exists ({info['residues']} residues)")
        return True
    
    try:
        print(f"Downloading {pdb_id}: {info['description']} ({info['residues']} residues)...")
        urllib.request.urlretrieve(url, output_file)
        print(f"✓ Downloaded {pdb_id}")
        return True
    except Exception as e:
        print(f"✗ Failed to download {pdb_id}: {e}")
        return False

def main():
    print("=" * 70)
    print("Downloading RINmaker Test Proteins")
    print("=" * 70)
    
    success_count = 0
    for pdb_id, info in TEST_PROTEINS.items():
        if download_pdb(pdb_id, info):
            success_count += 1
    
    print("\n" + "=" * 70)
    print(f"Downloaded {success_count}/{len(TEST_PROTEINS)} structures")
    print("=" * 70)

if __name__ == '__main__':
    main()
