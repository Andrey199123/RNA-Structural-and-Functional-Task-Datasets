#!/usr/bin/env python3
"""
Extract RNA base pairing and stacking interactions from PDB files.

Uses geometric criteria to detect:
1. Base pairs - based on hydrogen bonding distances between bases
2. Stacking interactions - based on vertical distance and ring parallelism

Output: NPZ files with N×N interaction matrices for each RNA.
"""

import argparse
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# Atoms used for base pairing detection (Watson-Crick edge)
WC_ATOMS = {
    'A': ['N1', 'N6'],
    'G': ['N1', 'O6', 'N2'],
    'C': ['N3', 'N4', 'O2'],
    'U': ['N3', 'O4', 'O2'],
    'DA': ['N1', 'N6'],
    'DG': ['N1', 'O6', 'N2'],
    'DC': ['N3', 'N4', 'O2'],
    'DT': ['N3', 'O4', 'O2'],
}

# Ring center atoms for stacking detection
RING_ATOMS = {
    'A': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],  # Purine
    'G': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
    'C': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],  # Pyrimidine
    'U': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
    'DA': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
    'DG': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
    'DC': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
    'DT': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
}


def get_base_type(resname: str) -> str:
    """Normalize residue name to single letter."""
    mapping = {'A': 'A', 'G': 'G', 'C': 'C', 'U': 'U',
               'DA': 'A', 'DG': 'G', 'DC': 'C', 'DT': 'U'}
    return mapping.get(resname, None)


def get_ring_center(residue) -> Optional[np.ndarray]:
    """Calculate center of base ring atoms."""
    resname = residue.get_resname()
    if resname not in RING_ATOMS:
        return None
    
    coords = []
    for atom_name in RING_ATOMS[resname]:
        if atom_name in residue:
            coords.append(residue[atom_name].get_coord())
    
    if len(coords) < 3:
        return None
    
    return np.mean(coords, axis=0)


def get_ring_normal(residue) -> Optional[np.ndarray]:
    """Calculate normal vector to base plane using ring atoms."""
    resname = residue.get_resname()
    if resname not in RING_ATOMS:
        return None
    
    coords = []
    for atom_name in RING_ATOMS[resname][:3]:  # Use first 3 atoms
        if atom_name in residue:
            coords.append(residue[atom_name].get_coord())
    
    if len(coords) < 3:
        return None
    
    # Calculate normal from cross product
    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[0]
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)
    if norm < 1e-6:
        return None
    return normal / norm


def detect_base_pair(res1, res2, max_dist: float = 3.5) -> bool:
    """
    Detect if two residues form a base pair.
    Uses hydrogen bonding distance criterion between WC edge atoms.
    """
    resname1, resname2 = res1.get_resname(), res2.get_resname()
    
    if resname1 not in WC_ATOMS or resname2 not in WC_ATOMS:
        return False
    
    # Check distances between potential H-bonding atoms
    for atom1_name in WC_ATOMS[resname1]:
        if atom1_name not in res1:
            continue
        coord1 = res1[atom1_name].get_coord()
        
        for atom2_name in WC_ATOMS[resname2]:
            if atom2_name not in res2:
                continue
            coord2 = res2[atom2_name].get_coord()
            
            dist = np.linalg.norm(coord1 - coord2)
            if dist < max_dist:
                return True
    
    return False


def detect_stacking(res1, res2, max_vertical_dist: float = 4.5, 
                    max_angle: float = 30.0) -> bool:
    """
    Detect if two bases are stacked.
    Criteria: ring centers ~3.4Å apart vertically, rings roughly parallel.
    """
    center1 = get_ring_center(res1)
    center2 = get_ring_center(res2)
    normal1 = get_ring_normal(res1)
    normal2 = get_ring_normal(res2)
    
    if center1 is None or center2 is None:
        return False
    if normal1 is None or normal2 is None:
        return False
    
    # Check vertical distance (projection onto normal)
    displacement = center2 - center1
    vertical_dist = abs(np.dot(displacement, normal1))
    
    if vertical_dist > max_vertical_dist or vertical_dist < 2.5:
        return False
    
    # Check horizontal displacement (should be small for stacking)
    horizontal_dist = np.sqrt(np.linalg.norm(displacement)**2 - vertical_dist**2)
    if horizontal_dist > 2.5:  # Centers should be roughly aligned
        return False
    
    # Check angle between ring planes
    cos_angle = abs(np.dot(normal1, normal2))
    angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
    
    return angle < max_angle


def extract_interactions(pdb_file: Path) -> Optional[Dict]:
    """Extract base pairing and stacking matrices from PDB."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('rna', pdb_file)
    except Exception:
        return None
    
    # Collect residues
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ':  # Standard residue
                    resname = residue.get_resname()
                    if resname in WC_ATOMS:
                        residues.append(residue)
        break  # First model only
    
    n = len(residues)
    if n < 2:
        return None
    
    # Initialize matrices
    base_pair_matrix = np.zeros((n, n), dtype=np.float32)
    stacking_matrix = np.zeros((n, n), dtype=np.float32)

    
    # Detect interactions
    for i in range(n):
        for j in range(i + 1, n):
            # Skip adjacent residues for base pairing (they can't pair)
            if abs(i - j) > 3:  # At least 3 residues apart for base pair
                if detect_base_pair(residues[i], residues[j]):
                    base_pair_matrix[i, j] = 1.0
                    base_pair_matrix[j, i] = 1.0
            
            # Stacking typically between adjacent or near-adjacent bases
            if abs(i - j) <= 2:
                if detect_stacking(residues[i], residues[j]):
                    stacking_matrix[i, j] = 1.0
                    stacking_matrix[j, i] = 1.0
    
    return {
        'base_pair_matrix': base_pair_matrix,
        'stacking_matrix': stacking_matrix,
        'sequence_length': n,
        'num_base_pairs': int(np.sum(base_pair_matrix) / 2),
        'num_stacking': int(np.sum(stacking_matrix) / 2),
    }


def process_single_file(args_tuple):
    """Process a single PDB file (for parallel execution)."""
    pdb_file, output_dir = args_tuple
    rna_id = pdb_file.stem
    output_file = output_dir / f"{rna_id}_interactions.npz"
    
    if output_file.exists():
        return (rna_id, True, 0, 0, 0, "cached")
    
    result = extract_interactions(pdb_file)
    
    if result is not None:
        np.savez_compressed(output_file, **result)
        return (rna_id, True, result['sequence_length'], 
                result['num_base_pairs'], result['num_stacking'], "ok")
    else:
        return (rna_id, False, 0, 0, 0, "failed")


def main():
    from concurrent.futures import ProcessPoolExecutor, as_completed
    
    parser = argparse.ArgumentParser(description="Extract RNA base interactions")
    parser.add_argument('--input-dir', type=Path,
                       default=Path('RNA-Protein/data/structural_features/rna_structures'),
                       help='Directory containing RNA PDB files')
    parser.add_argument('--output-dir', type=Path,
                       default=Path('RNA-Protein/data/structural_features/rna_interactions'),
                       help='Output directory for NPZ files')
    parser.add_argument('--workers', type=int, default=8,
                       help='Number of parallel workers')
    args = parser.parse_args()
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    pdb_files = list(args.input_dir.glob('*.pdb'))
    print(f"Found {len(pdb_files)} PDB files")
    print(f"Using {args.workers} workers")
    
    work_items = [(pdb, args.output_dir) for pdb in pdb_files]
    
    success = 0
    failed = 0
    total_bp = 0
    total_stack = 0
    processed = 0
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_single_file, item): item for item in work_items}
        
        for future in as_completed(futures):
            processed += 1
            rna_id, ok, seq_len, num_bp, num_stack, status = future.result()
            
            if ok:
                success += 1
                total_bp += num_bp
                total_stack += num_stack
                if processed % 100 == 0 or processed == len(pdb_files):
                    print(f"[{processed}/{len(pdb_files)}] {rna_id}: {seq_len}nt, {num_bp} bp, {num_stack} stacks")
            else:
                failed += 1
                print(f"[{processed}/{len(pdb_files)}] Failed: {rna_id}")
    
    print(f"\nDone! Success: {success}, Failed: {failed}")
    print(f"Total base pairs detected: {total_bp}")
    print(f"Total stacking interactions: {total_stack}")


if __name__ == "__main__":
    main()
