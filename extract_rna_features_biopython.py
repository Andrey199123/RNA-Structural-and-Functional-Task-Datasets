#!/usr/bin/env python3
"""
Extract RNA structural features using Biopython (DSSR alternative)
"""

from Bio.PDB import PDBParser, NeighborSearch
import numpy as np
from pathlib import Path
import json

def extract_rna_features(pdb_file):
    """
    Extract structural features from RNA PDB file using Biopython
    
    Args:
        pdb_file: Path to PDB file
        
    Returns:
        Dictionary of features
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('RNA', pdb_file)
    
    features = {}
    
    # Get all residues and atoms
    residues = list(structure.get_residues())
    atoms = list(structure.get_atoms())
    
    # Basic counts
    features['num_residues'] = len(residues)
    features['num_chains'] = len(list(structure.get_chains()))
    features['num_atoms'] = len(atoms)
    
    # Nucleotide composition
    residue_names = [res.resname.strip() for res in residues]
    features['num_A'] = residue_names.count('A')
    features['num_U'] = residue_names.count('U') + residue_names.count('T')
    features['num_G'] = residue_names.count('G')
    features['num_C'] = residue_names.count('C')
    
    # Sequence length
    features['sequence_length'] = len(residue_names)
    
    # GC content
    gc_count = features['num_G'] + features['num_C']
    features['gc_content'] = gc_count / len(residue_names) if residue_names else 0
    
    # Geometric features
    coords = np.array([atom.coord for atom in atoms])
    
    if len(coords) > 0:
        # Center of mass
        center = coords.mean(axis=0)
        features['center_x'] = float(center[0])
        features['center_y'] = float(center[1])
        features['center_z'] = float(center[2])
        
        # Radius of gyration
        features['radius_gyration'] = float(np.sqrt(((coords - center)**2).sum(axis=1).mean()))
        
        # Bounding box dimensions
        features['bbox_x'] = float(coords[:, 0].max() - coords[:, 0].min())
        features['bbox_y'] = float(coords[:, 1].max() - coords[:, 1].min())
        features['bbox_z'] = float(coords[:, 2].max() - coords[:, 2].min())
        features['bbox_volume'] = features['bbox_x'] * features['bbox_y'] * features['bbox_z']
    
    # Estimate base pairs (simplified - distance-based)
    features['estimated_base_pairs'] = estimate_base_pairs(residues)
    
    # Secondary structure estimate
    features['estimated_stems'] = estimate_stems(residues)
    
    return features


def estimate_base_pairs(residues):
    """
    Estimate number of base pairs based on C1' atom distances
    Watson-Crick pairs have C1'-C1' distance ~10.5 Å
    """
    c1_atoms = []
    for res in residues:
        if "C1'" in res:
            c1_atoms.append(res["C1'"])
    
    if len(c1_atoms) < 2:
        return 0
    
    # Use NeighborSearch to find close pairs
    ns = NeighborSearch(c1_atoms)
    pairs = 0
    
    for atom in c1_atoms:
        # Find atoms within 9-12 Å (typical base pair distance)
        neighbors = ns.search(atom.coord, 12.0, 'A')
        for neighbor in neighbors:
            if neighbor != atom:
                dist = atom - neighbor
                if 9.0 < dist < 12.0:  # Typical base pair distance
                    pairs += 1
    
    # Divide by 2 since we count each pair twice
    return pairs // 2


def estimate_stems(residues):
    """
    Rough estimate of stem regions based on consecutive base pairs
    """
    # This is a simplified heuristic
    # Real stem detection would need proper base pair identification
    num_residues = len(residues)
    
    if num_residues < 10:
        return 0
    elif num_residues < 30:
        return 1
    elif num_residues < 60:
        return 2
    else:
        return 3


def test_extraction():
    """Test feature extraction on example RNA structure"""
    
    print("=" * 60)
    print("Testing RNA Feature Extraction (Biopython)")
    print("=" * 60)
    print()
    
    # Test files
    test_files = [
        "test_output/example_structure/model_1_unrelaxed.pdb",
        "test_output/single_structure/model_1_unrelaxed.pdb"
    ]
    
    for pdb_file in test_files:
        pdb_path = Path(pdb_file)
        
        if not pdb_path.exists():
            print(f"✗ File not found: {pdb_file}")
            continue
        
        print(f"Analyzing: {pdb_path.name}")
        print("-" * 60)
        
        try:
            features = extract_rna_features(pdb_path)
            
            # Print features
            print(f"  Residues: {features['num_residues']}")
            print(f"  Chains: {features['num_chains']}")
            print(f"  Atoms: {features['num_atoms']}")
            print(f"  Composition: A={features['num_A']}, U={features['num_U']}, "
                  f"G={features['num_G']}, C={features['num_C']}")
            print(f"  GC content: {features['gc_content']:.2%}")
            print(f"  Radius of gyration: {features['radius_gyration']:.2f} Å")
            print(f"  Estimated base pairs: {features['estimated_base_pairs']}")
            print(f"  Estimated stems: {features['estimated_stems']}")
            
            # Save to JSON
            output_json = pdb_path.parent / f"{pdb_path.stem}_features.json"
            with open(output_json, 'w') as f:
                json.dump(features, f, indent=2)
            
            print(f"  ✓ Features saved to: {output_json}")
            print()
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
            print()
    
    print("=" * 60)
    print("✓ Feature extraction working!")
    print()
    print("This is a simplified alternative to DSSR.")
    print("Features are less detailed but sufficient for ML model.")
    print("=" * 60)


if __name__ == "__main__":
    test_extraction()
