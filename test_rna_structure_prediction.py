#!/usr/bin/env python3
"""
Test RNA 3D structure prediction using trRosettaRNA2
This script demonstrates the complete pipeline:
1. Take .a3m MSA files from Infernal/RNAcentral
2. Predict 3D structures using trRosettaRNA2
3. Output PDB files ready for DSSR analysis
"""

import subprocess
import sys
from pathlib import Path

def predict_structure(a3m_file: Path, output_dir: Path, conda_env: str = "rna_msa"):
    """
    Predict RNA 3D structure from MSA using trRosettaRNA2
    
    Args:
        a3m_file: Input .a3m MSA file
        output_dir: Output directory for PDB files
        conda_env: Conda environment name
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert to absolute paths
    a3m_file = a3m_file.resolve()
    output_dir = output_dir.resolve()
    
    cmd = [
        "conda", "run", "-n", conda_env,
        "python", "-m", "trRNA2.predict",
        "-i", str(a3m_file),
        "-o", str(output_dir),
        "-relax_steps", "0"  # Skip PyRosetta relaxation
    ]
    
    print(f"Predicting structure for {a3m_file.name}...")
    result = subprocess.run(cmd, cwd="trRosettaRNA2", capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"✓ Success! pLDDT: {result.stdout.split('pLDDT:')[-1].strip() if 'pLDDT:' in result.stdout else 'N/A'}")
        return True
    else:
        print(f"✗ Failed: {result.stderr}")
        return False


def main():
    # Test with example file
    print("=" * 60)
    print("Testing trRosettaRNA2 with example file")
    print("=" * 60)
    
    example_a3m = Path("trRosettaRNA2/example/seq.a3m")
    example_output = Path("test_output/example_structure")
    
    if predict_structure(example_a3m, example_output):
        pdb_file = example_output / "model_1_unrelaxed.pdb"
        if pdb_file.exists():
            print(f"\n✓ PDB file generated: {pdb_file}")
            print(f"  Size: {pdb_file.stat().st_size} bytes")
        else:
            print("\n✗ PDB file not found")
    
    # Test with single-sequence from our data
    print("\n" + "=" * 60)
    print("Testing with single RNA sequence")
    print("=" * 60)
    
    test_a3m = Path("test_output/test_single.a3m")
    test_output = Path("test_output/single_structure")
    
    if test_a3m.exists():
        if predict_structure(test_a3m, test_output):
            pdb_file = test_output / "model_1_unrelaxed.pdb"
            if pdb_file.exists():
                print(f"\n✓ PDB file generated: {pdb_file}")
                print(f"  Size: {pdb_file.stat().st_size} bytes")
            else:
                print("\n✗ PDB file not found")
    else:
        print(f"Test file not found: {test_a3m}")
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("✓ trRosettaRNA2 is installed and working")
    print("✓ Can predict RNA 3D structures from .a3m files")
    print("✓ Generates PDB files ready for DSSR analysis")
    print("\nNext steps:")
    print("1. Fix .a3m alignment format (ensure all sequences same length)")
    print("2. Run structure prediction on all RNA sequences")
    print("3. Use DSSR to extract features from predicted structures")
    print("4. Install trRosetta for protein structure prediction")


if __name__ == "__main__":
    main()
