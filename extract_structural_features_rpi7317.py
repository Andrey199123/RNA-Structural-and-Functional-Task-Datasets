#!/usr/bin/env python3
"""
Extract structural features for RPI7317 dataset

This script processes the RPI7317 dataset to generate structural features:
1. Generate protein structures using ESMFold API
2. Extract protein features using RINmaker
3. Generate RNA structures using trRosettaRNA2 (with MSAs)
4. Extract RNA features using Biopython

Usage:
    python extract_structural_features_rpi7317.py --mode test  # Process 10 samples
    python extract_structural_features_rpi7317.py --mode full  # Process all samples
"""

import argparse
import json
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Tuple
import requests
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
import numpy as np


class StructuralFeatureExtractor:
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.protein_pdb_dir = output_dir / "protein_structures"
        self.rna_pdb_dir = output_dir / "rna_structures"
        self.protein_features_dir = output_dir / "protein_features"
        self.rna_features_dir = output_dir / "rna_features"
        
        # Create directories
        for d in [self.protein_pdb_dir, self.rna_pdb_dir, 
                  self.protein_features_dir, self.rna_features_dir]:
            d.mkdir(parents=True, exist_ok=True)
        
        self.rinmaker_path = Path("RINmaker/RINmaker/build/app/RINmaker")
        
    def load_sequences(self, fasta_path: Path) -> Dict[str, str]:
        """Load sequences from FASTA file"""
        sequences = {}
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Extract ID (handle both protein and RNA formats)
            seq_id = record.id.split("|")[1] if "|" in record.id else record.id
            sequences[seq_id] = str(record.seq)
        return sequences
    
    def predict_protein_structure_esmfold(self, protein_id: str, sequence: str) -> Path:
        """Predict protein structure using ESMFold API"""
        output_file = self.protein_pdb_dir / f"{protein_id}.pdb"
        
        if output_file.exists():
            print(f"  Protein structure already exists: {protein_id}")
            return output_file
        
        print(f"  Predicting protein structure: {protein_id} ({len(sequence)} aa)")
        
        # ESMFold API has a length limit of ~400 amino acids
        if len(sequence) > 400:
            print(f"    ✗ Sequence too long for ESMFold API (max ~400 aa)")
            return None
        
        try:
            response = requests.post(
                'https://api.esmatlas.com/foldSequence/v1/pdb/',
                data=sequence,
                headers={'Content-Type': 'text/plain'},
                timeout=60
            )
            
            if response.status_code == 200:
                output_file.write_text(response.text)
                print(f"    ✓ Saved to {output_file.name}")
                return output_file
            else:
                print(f"    ✗ API error: {response.status_code}")
                return None
                
        except Exception as e:
            print(f"    ✗ Error: {e}")
            return None
    
    def extract_protein_features_rinmaker(self, protein_id: str, pdb_file: Path) -> Dict:
        """Extract protein features using RINmaker"""
        output_file = self.protein_features_dir / f"{protein_id}_features.json"
        
        if output_file.exists():
            print(f"  Protein features already exist: {protein_id}")
            return json.loads(output_file.read_text())
        
        print(f"  Extracting protein features: {protein_id}")
        
        if not self.rinmaker_path.exists():
            print(f"    ✗ RINmaker not found at {self.rinmaker_path}")
            return None
        
        try:
            # Run RINmaker
            output_base = self.protein_features_dir / f"{protein_id}"
            graphml_file = self.protein_features_dir / f"{protein_id}.graphml"
            result = subprocess.run(
                [str(self.rinmaker_path), "--input", str(pdb_file), "--output", str(output_base), "rin"],
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode != 0:
                print(f"    ✗ RINmaker failed: {result.stderr}")
                return None
            
            # Check if output file exists
            if not graphml_file.exists():
                print(f"    ✗ Output file not created: {graphml_file}")
                return None
            
            # Parse GraphML and extract basic features
            # For now, just count nodes and edges
            graphml_content = graphml_file.read_text()
            num_nodes = graphml_content.count('<node id=')
            num_edges = graphml_content.count('<edge ')
            
            features = {
                "protein_id": protein_id,
                "num_residues": num_nodes,
                "num_interactions": num_edges,
                "interaction_density": num_edges / (num_nodes * (num_nodes - 1)) if num_nodes > 1 else 0,
                "graphml_file": str(graphml_file)
            }
            
            output_file.write_text(json.dumps(features, indent=2))
            print(f"    ✓ Features: {num_nodes} residues, {num_edges} interactions")
            return features
            
        except Exception as e:
            print(f"    ✗ Error: {e}")
            return None
    
    def extract_rna_features_biopython(self, rna_id: str, pdb_file: Path) -> Dict:
        """Extract RNA features using Biopython"""
        output_file = self.rna_features_dir / f"{rna_id}_features.json"
        
        if output_file.exists():
            print(f"  RNA features already exist: {rna_id}")
            return json.loads(output_file.read_text())
        
        print(f"  Extracting RNA features: {rna_id}")
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(rna_id, pdb_file)
            
            # Extract coordinates
            coords = []
            nucleotides = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[0] == ' ':  # Standard residue
                            nucleotides.append(residue.get_resname())
                            for atom in residue:
                                coords.append(atom.get_coord())
            
            coords = np.array(coords)
            
            # Calculate features
            center = coords.mean(axis=0)
            distances = np.linalg.norm(coords - center, axis=1)
            radius_of_gyration = np.sqrt(np.mean(distances**2))
            
            bbox_min = coords.min(axis=0)
            bbox_max = coords.max(axis=0)
            bbox_size = bbox_max - bbox_min
            
            # Nucleotide composition
            nuc_counts = {
                'A': nucleotides.count('A') + nucleotides.count('DA'),
                'U': nucleotides.count('U') + nucleotides.count('DT'),
                'G': nucleotides.count('G') + nucleotides.count('DG'),
                'C': nucleotides.count('C') + nucleotides.count('DC')
            }
            total_nucs = sum(nuc_counts.values())
            gc_content = (nuc_counts['G'] + nuc_counts['C']) / total_nucs if total_nucs > 0 else 0
            
            features = {
                "rna_id": rna_id,
                "num_nucleotides": len(nucleotides),
                "num_atoms": len(coords),
                "radius_of_gyration": float(radius_of_gyration),
                "bounding_box": {
                    "x": float(bbox_size[0]),
                    "y": float(bbox_size[1]),
                    "z": float(bbox_size[2]),
                    "volume": float(np.prod(bbox_size))
                },
                "nucleotide_composition": nuc_counts,
                "gc_content": float(gc_content)
            }
            
            output_file.write_text(json.dumps(features, indent=2))
            print(f"    ✓ Features: {len(nucleotides)} nucleotides, GC={gc_content:.2f}")
            return features
            
        except Exception as e:
            print(f"    ✗ Error: {e}")
            return None


def main():
    parser = argparse.ArgumentParser(description="Extract structural features for RPI7317")
    parser.add_argument("--mode", choices=["test", "full"], default="test",
                       help="test: process 10 samples, full: process all")
    parser.add_argument("--output-dir", type=Path, default=Path("RNA-Protein/data/structural_features"),
                       help="Output directory for features")
    args = parser.parse_args()
    
    print("=" * 80)
    print("RPI7317 Structural Feature Extraction Pipeline")
    print("=" * 80)
    print()
    
    # Initialize extractor
    extractor = StructuralFeatureExtractor(args.output_dir)
    
    # Load sequences
    print("Loading sequences...")
    protein_fasta = Path("RNA-Protein/data/RPI7317/RPI7317_protein_sequence.fasta")
    rna_fasta = Path("RNA-Protein/data/RPI7317/RPI7317_IncRNA_sequence.fasta")
    
    proteins = extractor.load_sequences(protein_fasta)
    rnas = extractor.load_sequences(rna_fasta)
    
    print(f"  Loaded {len(proteins)} proteins")
    print(f"  Loaded {len(rnas)} RNAs")
    print()
    
    # Limit to test samples if in test mode
    if args.mode == "test":
        proteins = dict(list(proteins.items())[:10])
        rnas = dict(list(rnas.items())[:10])
        print(f"TEST MODE: Processing {len(proteins)} proteins and {len(rnas)} RNAs")
        print()
    
    # Process proteins
    print("=" * 80)
    print("STEP 1: Generating Protein Structures (ESMFold API)")
    print("=" * 80)
    print()
    
    protein_stats = {"success": 0, "failed": 0}
    for i, (prot_id, sequence) in enumerate(proteins.items(), 1):
        print(f"[{i}/{len(proteins)}] Processing protein: {prot_id}")
        pdb_file = extractor.predict_protein_structure_esmfold(prot_id, sequence)
        if pdb_file:
            protein_stats["success"] += 1
        else:
            protein_stats["failed"] += 1
        time.sleep(1)  # Rate limiting
        print()
    
    print(f"Protein structure generation: {protein_stats['success']} success, {protein_stats['failed']} failed")
    print()
    
    # Extract protein features
    print("=" * 80)
    print("STEP 2: Extracting Protein Features (RINmaker)")
    print("=" * 80)
    print()
    
    protein_feature_stats = {"success": 0, "failed": 0}
    for i, prot_id in enumerate(proteins.keys(), 1):
        pdb_file = extractor.protein_pdb_dir / f"{prot_id}.pdb"
        if not pdb_file.exists():
            print(f"[{i}/{len(proteins)}] Skipping {prot_id} (no PDB file)")
            protein_feature_stats["failed"] += 1
            continue
        
        print(f"[{i}/{len(proteins)}] Processing protein: {prot_id}")
        features = extractor.extract_protein_features_rinmaker(prot_id, pdb_file)
        if features:
            protein_feature_stats["success"] += 1
        else:
            protein_feature_stats["failed"] += 1
        print()
    
    print(f"Protein feature extraction: {protein_feature_stats['success']} success, {protein_feature_stats['failed']} failed")
    print()
    
    # Note about RNA structures
    print("=" * 80)
    print("STEP 3: RNA Structure Generation (Manual)")
    print("=" * 80)
    print()
    print("RNA structure generation requires:")
    print("  1. MSA generation using your hybrid pipeline (Infernal+Rfam + BLAST+RNAcentral)")
    print("  2. trRosettaRNA2 prediction (requires conda environment 'rna_msa')")
    print()
    print("This step is NOT automated in this script due to:")
    print("  - Long runtime (hours for 1,874 RNAs)")
    print("  - Requires conda environment activation")
    print("  - Requires large databases (RNAcentral, Rfam)")
    print()
    print("To generate RNA structures manually:")
    print("  1. Activate conda environment: conda activate rna_msa")
    print("  2. Run your hybrid MSA pipeline on RPI7317 RNAs")
    print("  3. Run trRosettaRNA2 on each MSA to generate PDB files")
    print("  4. Place PDB files in: RNA-Protein/data/structural_features/rna_structures/")
    print()
    print("Once RNA PDB files are available, run this script again to extract features.")
    print()
    
    # Extract RNA features (if PDB files exist)
    print("=" * 80)
    print("STEP 4: Extracting RNA Features (Biopython)")
    print("=" * 80)
    print()
    
    rna_feature_stats = {"success": 0, "failed": 0, "no_pdb": 0}
    for i, rna_id in enumerate(rnas.keys(), 1):
        pdb_file = extractor.rna_pdb_dir / f"{rna_id}.pdb"
        if not pdb_file.exists():
            print(f"[{i}/{len(rnas)}] Skipping {rna_id} (no PDB file)")
            rna_feature_stats["no_pdb"] += 1
            continue
        
        print(f"[{i}/{len(rnas)}] Processing RNA: {rna_id}")
        features = extractor.extract_rna_features_biopython(rna_id, pdb_file)
        if features:
            rna_feature_stats["success"] += 1
        else:
            rna_feature_stats["failed"] += 1
        print()
    
    print(f"RNA feature extraction: {rna_feature_stats['success']} success, {rna_feature_stats['failed']} failed, {rna_feature_stats['no_pdb']} no PDB")
    print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print(f"Proteins processed: {len(proteins)}")
    print(f"  - Structures generated: {protein_stats['success']}")
    print(f"  - Features extracted: {protein_feature_stats['success']}")
    print()
    print(f"RNAs processed: {len(rnas)}")
    print(f"  - Structures available: {rna_feature_stats['success'] + rna_feature_stats['failed']}")
    print(f"  - Features extracted: {rna_feature_stats['success']}")
    print()
    print(f"Output directory: {args.output_dir}")
    print()


if __name__ == "__main__":
    main()
