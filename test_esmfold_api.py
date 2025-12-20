#!/usr/bin/env python3
"""
Test ESMFold API for protein structure prediction
No local installation needed - uses Meta's free API
"""

import requests
import time

# Test protein sequence (small for quick test)
test_sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"

def predict_structure_api(sequence, output_file=None):
    """
    Predict protein structure using ESMFold API
    
    Args:
        sequence: Protein amino acid sequence
        output_file: Optional path to save PDB file
        
    Returns:
        PDB string
    """
    print(f"Predicting structure for {len(sequence)} amino acid sequence...")
    print("Using ESMFold API (https://api.esmatlas.com)")
    
    try:
        response = requests.post(
            'https://api.esmatlas.com/foldSequence/v1/pdb/',
            data=sequence,
            headers={'Content-Type': 'text/plain'},
            timeout=60
        )
        
        if response.status_code == 200:
            pdb_string = response.text
            
            if output_file:
                with open(output_file, 'w') as f:
                    f.write(pdb_string)
                print(f"✓ Saved to: {output_file}")
            
            return pdb_string
        else:
            print(f"✗ API error: {response.status_code}")
            print(response.text)
            return None
            
    except Exception as e:
        print(f"✗ Request failed: {e}")
        return None


def main():
    import os
    os.makedirs("test_output", exist_ok=True)
    
    print("=" * 60)
    print("Testing ESMFold API")
    print("=" * 60)
    print()
    
    start_time = time.time()
    
    pdb_string = predict_structure_api(
        test_sequence,
        output_file="test_output/esmfold_api_test.pdb"
    )
    
    elapsed = time.time() - start_time
    
    if pdb_string:
        print()
        print("✓ Prediction successful!")
        print(f"  Time: {elapsed:.1f} seconds")
        print(f"  PDB size: {len(pdb_string)} characters")
        print()
        print("PDB preview (first 500 chars):")
        print(pdb_string[:500])
        print()
        print("=" * 60)
        print("✓ ESMFold API is working!")
        print()
        print("Usage for your project:")
        print("  1. Get protein sequences from your dataset")
        print("  2. Call predict_structure_api() for each sequence")
        print("  3. Save PDB files")
        print("  4. Use RINmaker to extract features from PDB files")
        print("=" * 60)
    else:
        print()
        print("=" * 60)
        print("✗ ESMFold API test failed")
        print("Check your internet connection")
        print("=" * 60)


if __name__ == "__main__":
    main()
