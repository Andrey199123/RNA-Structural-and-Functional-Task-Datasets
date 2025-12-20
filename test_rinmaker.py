#!/usr/bin/env python3
"""
Test RINmaker for protein feature extraction
"""

import subprocess
import os
from pathlib import Path

def test_rinmaker():
    """Test RINmaker on a protein PDB file"""
    
    print("=" * 60)
    print("Testing RINmaker")
    print("=" * 60)
    print()
    
    # Path to RINmaker binary
    rinmaker_bin = Path("RINmaker/RINmaker/build/app/RINmaker")
    
    if not rinmaker_bin.exists():
        print("✗ RINmaker binary not found")
        print(f"  Expected at: {rinmaker_bin}")
        return False
    
    # Test with ESMFold predicted structure
    test_pdb = Path("test_output/esmfold_api_test.pdb")
    
    if not test_pdb.exists():
        print(f"✗ Test PDB file not found: {test_pdb}")
        print("  Run test_esmfold_api.py first to generate a test structure")
        return False
    
    output_file = Path("test_output/rinmaker_test.graphml")
    
    print(f"Input PDB: {test_pdb}")
    print(f"Output: {output_file}")
    print()
    print("Running RINmaker...")
    
    try:
        result = subprocess.run(
            [str(rinmaker_bin), "-i", str(test_pdb), "-o", str(output_file), "rin"],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        # RINmaker adds .graphml extension, so check for that
        actual_output = Path(str(output_file) + ".graphml")
        
        if actual_output.exists():
            size = actual_output.stat().st_size
            print(f"✓ RINmaker completed successfully!")
            print(f"  Output file: {actual_output}")
            print(f"  Size: {size:,} bytes")
            
            # Count nodes and edges
            with open(actual_output) as f:
                content = f.read()
                nodes = content.count("<node id=")
                edges = content.count("<edge ")
            
            print(f"  Nodes (residues): {nodes}")
            print(f"  Edges (interactions): {edges}")
            
            return True
        else:
            print(f"✗ Output file not created")
            print(f"  Stdout: {result.stdout}")
            print(f"  Stderr: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ RINmaker timed out")
        return False
    except Exception as e:
        print(f"✗ Error running RINmaker: {e}")
        return False


def main():
    os.makedirs("test_output", exist_ok=True)
    
    success = test_rinmaker()
    
    print()
    print("=" * 60)
    if success:
        print("✓ RINmaker is working!")
        print()
        print("Usage for your project:")
        print("  1. Get protein PDB files (from ESMFold API)")
        print("  2. Run RINmaker to extract interaction networks")
        print("  3. Parse GraphML output to get features")
        print("  4. Use features in ML model")
        print()
        print("Command:")
        print("  ./RINmaker/RINmaker/build/app/RINmaker \\")
        print("    -i protein.pdb -o output.graphml rin")
    else:
        print("✗ RINmaker test failed")
    print("=" * 60)


if __name__ == "__main__":
    main()
