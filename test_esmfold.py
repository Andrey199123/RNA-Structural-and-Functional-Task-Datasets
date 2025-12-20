#!/usr/bin/env python3
"""
Test ESMFold for protein structure prediction
"""

import torch
import esm

# Test sequence (small protein for quick test)
test_sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"

def test_esmfold():
    """Test ESMFold structure prediction"""
    
    print("=" * 60)
    print("Testing ESMFold Installation")
    print("=" * 60)
    
    # Check if CUDA is available
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\nDevice: {device}")
    
    # Load ESMFold model
    print("\nLoading ESMFold model...")
    print("(This will download ~2GB model weights on first run)")
    
    try:
        model = esm.pretrained.esmfold_v1()
        model = model.eval()
        if device == "cuda":
            model = model.cuda()
        print("✓ Model loaded successfully")
    except Exception as e:
        print(f"✗ Failed to load model: {e}")
        return False
    
    # Test prediction
    print(f"\nTesting prediction on {len(test_sequence)} amino acid sequence...")
    
    try:
        with torch.no_grad():
            output = model.infer_pdb(test_sequence)
        
        print("✓ Prediction successful!")
        print(f"\nPDB output preview (first 500 chars):")
        print(output[:500])
        
        # Save test output
        with open("test_output/esmfold_test.pdb", "w") as f:
            f.write(output)
        print(f"\n✓ Saved to: test_output/esmfold_test.pdb")
        
        return True
        
    except Exception as e:
        print(f"✗ Prediction failed: {e}")
        return False

if __name__ == "__main__":
    import os
    os.makedirs("test_output", exist_ok=True)
    
    success = test_esmfold()
    
    print("\n" + "=" * 60)
    if success:
        print("✓ ESMFold is working!")
        print("\nUsage:")
        print("  import esm")
        print("  model = esm.pretrained.esmfold_v1()")
        print("  pdb_string = model.infer_pdb(protein_sequence)")
    else:
        print("✗ ESMFold test failed")
    print("=" * 60)
