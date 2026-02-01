#!/usr/bin/env python3
"""
Extract 2D structural features from protein PDB files.

Generates distance maps and contact maps that can be used as
image-like inputs for the model (similar to trRosetta approach).

Features extracted:
1. Distance matrix (L x L) - Cα atom distances between residues
2. Contact map (L x L) - Binary contacts at 8Å threshold
3. Normalized distance map (L x L) - Distances scaled to [0, 1]

Output: NPZ files containing these matrices for each protein.
"""

import argparse
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser
from typing import Dict, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')


def extract_ca_coords(pdb_file: Path) -> Tuple[np.ndarray, list]:
    """
    Extract Cα atom coordinates from protein PDB file.
    Cα is the standard reference atom for proteins.
    
    Returns:
        coords: (L, 3) array of Cα coordinates
        residue_ids: list of residue identifiers
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    coords = []
    residue_ids = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip hetero atoms and water
                if residue.get_id()[0] != ' ':
                    continue
                
                # Get Cα atom
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
                    residue_ids.append(f"{chain.id}_{residue.get_id()[1]}")
        break  # Only first model
    
    return np.array(coords), residue_ids


def compute_distance_matrix(coords: np.ndarray) -> np.ndarray:
    """Compute pairwise distance matrix from coordinates."""
    n = len(coords)
    dist_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            dist_matrix[i, j] = np.linalg.norm(coords[i] - coords[j])
    
    return dist_matrix


def compute_contact_map(dist_matrix: np.ndarray, threshold: float = 8.0) -> np.ndarray:
    """Compute binary contact map from distance matrix."""
    return (dist_matrix < threshold).astype(np.float32)


def normalize_distance_matrix(dist_matrix: np.ndarray, max_dist: float = 40.0) -> np.ndarray:
    """Normalize distance matrix to [0, 1] range."""
    normalized = np.clip(dist_matrix / max_dist, 0, 1)
    return normalized.astype(np.float32)


def extract_features(pdb_file: Path) -> Optional[Dict[str, np.ndarray]]:
    """Extract all 2D features from a single PDB file."""
    try:
        coords, residue_ids = extract_ca_coords(pdb_file)
        
        if len(coords) < 2:
            return None
        
        # Compute features
        dist_matrix = compute_distance_matrix(coords)
        contact_map = compute_contact_map(dist_matrix, threshold=8.0)
        norm_dist = normalize_distance_matrix(dist_matrix, max_dist=40.0)
        
        return {
            'distance_matrix': dist_matrix.astype(np.float32),
            'contact_map': contact_map,
            'normalized_distance': norm_dist,
            'sequence_length': len(coords),
            'residue_ids': np.array(residue_ids, dtype=object)
        }
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Extract 2D structural features from protein PDBs")
    parser.add_argument('--input-dir', type=Path, 
                       default=Path('RNA-Protein/data/structural_features/protein_structures'),
                       help='Directory containing protein PDB files')
    parser.add_argument('--output-dir', type=Path,
                       default=Path('RNA-Protein/data/structural_features/protein_2d_features'),
                       help='Output directory for NPZ files')
    args = parser.parse_args()
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all PDB files
    pdb_files = list(args.input_dir.glob('*.pdb'))
    print(f"Found {len(pdb_files)} PDB files")
    
    success = 0
    failed = 0
    
    for i, pdb_file in enumerate(pdb_files, 1):
        protein_id = pdb_file.stem
        output_file = args.output_dir / f"{protein_id}_2d.npz"
        
        # Skip if already processed
        if output_file.exists():
            success += 1
            continue
        
        features = extract_features(pdb_file)
        
        if features is not None:
            np.savez_compressed(output_file, **features)
            success += 1
            print(f"[{i}/{len(pdb_files)}] {protein_id}: {features['sequence_length']} residues")
        else:
            failed += 1
            print(f"[{i}/{len(pdb_files)}] Failed: {protein_id}")
    
    print(f"\nDone! Success: {success}, Failed: {failed}")
    print(f"Output saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
