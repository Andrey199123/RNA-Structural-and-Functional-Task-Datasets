#!/usr/bin/env python3
"""
Test script for RNA MSA pipeline

This script runs a minimal test to verify the pipeline is working correctly.
"""

import sys
from pathlib import Path
import tempfile
import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from rna_msa.pipeline import RNASequence, RNAMSAPipeline
from rna_msa.msa_utils import check_dependencies
from rna_msa.feature_extraction import extract_2d_features, structure_to_features
from rna_msa.single_sequence import predict_structure_spotRNA


def test_dependencies():
    """Test 1: Check dependencies"""
    print("Test 1: Checking dependencies...")
    deps = check_dependencies()
    
    for tool, available in deps.items():
        status = "✓" if available else "✗"
        print(f"  {status} {tool}")
    
    # At minimum, we need numpy/scipy
    if not any(deps.values()):
        print("  ⚠ Warning: No external tools found, will use fallback methods")
    
    print("  ✓ Test 1 passed\n")
    return True


def test_rna_sequence():
    """Test 2: RNA sequence loading"""
    print("Test 2: Testing RNA sequence loading...")
    
    # Create test FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">test_rna_1\n")
        f.write("AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC\n")
        f.write(">test_rna_2\n")
        f.write("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC\n")
        temp_fasta = f.name
    
    try:
        sequences = RNASequence.from_fasta(Path(temp_fasta))
        
        assert len(sequences) == 2, f"Expected 2 sequences, got {len(sequences)}"
        assert sequences[0].id == "test_rna_1"
        assert sequences[0].length == 32
        assert sequences[1].id == "test_rna_2"
        
        print(f"  ✓ Loaded {len(sequences)} sequences")
        print(f"  ✓ Test 2 passed\n")
        return True
        
    finally:
        Path(temp_fasta).unlink()


def test_single_sequence_prediction():
    """Test 3: Single-sequence structure prediction"""
    print("Test 3: Testing single-sequence structure prediction...")
    
    # Create test RNA
    rna = RNASequence(
        id="test_rna",
        sequence="GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
        length=32
    )
    
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            structure = predict_structure_spotRNA(rna, Path(tmpdir))
            
            assert 'sequence' in structure
            assert 'structure' in structure
            assert len(structure['structure']) == rna.length
            
            print(f"  ✓ Predicted structure: {structure['structure'][:20]}...")
            print(f"  ✓ Method: {structure.get('method', 'unknown')}")
            print(f"  ✓ Test 3 passed\n")
            return True
            
        except Exception as e:
            print(f"  ✗ Test 3 failed: {e}\n")
            return False


def test_feature_extraction():
    """Test 4: Feature extraction from structure"""
    print("Test 4: Testing feature extraction...")
    
    # Create test structure
    structure_dict = {
        'sequence': 'GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC',
        'structure': '((((((((................))))))))',
        'mfe': -10.5,
        'accessibility': np.array([0.3] * 8 + [1.0] * 16 + [0.3] * 8)
    }
    
    try:
        features = structure_to_features(structure_dict)
        
        assert 'conservation' in features
        assert 'coevolution' in features
        assert 'accessibility' in features
        
        length = len(structure_dict['sequence'])
        assert features['conservation'].shape == (length,)
        assert features['coevolution'].shape == (length, length)
        assert features['accessibility'].shape == (length,)
        
        print(f"  ✓ Conservation shape: {features['conservation'].shape}")
        print(f"  ✓ Coevolution shape: {features['coevolution'].shape}")
        print(f"  ✓ Accessibility shape: {features['accessibility'].shape}")
        print(f"  ✓ Test 4 passed\n")
        return True
        
    except Exception as e:
        print(f"  ✗ Test 4 failed: {e}\n")
        return False


def test_pipeline_integration():
    """Test 5: Full pipeline integration"""
    print("Test 5: Testing full pipeline integration...")
    
    # Create test RNA
    rna = RNASequence(
        id="test_rna",
        sequence="AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC",
        length=48
    )
    
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            # Initialize pipeline
            pipeline = RNAMSAPipeline(
                output_dir=Path(tmpdir),
                n_cpu=1
            )
            
            # Process RNA (will use single-sequence fallback)
            result = pipeline.process_rna(rna, force_strategy=False)
            
            assert result['success'], f"Pipeline failed: {result.get('error', 'unknown')}"
            assert result['rna_id'] == rna.id
            assert 'strategy' in result
            
            print(f"  ✓ Pipeline completed successfully")
            print(f"  ✓ Strategy used: {result['strategy']}")
            print(f"  ✓ Test 5 passed\n")
            return True
            
        except Exception as e:
            print(f"  ✗ Test 5 failed: {e}\n")
            return False


def test_feature_saving():
    """Test 6: Feature saving and loading"""
    print("Test 6: Testing feature saving and loading...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            # Create test features
            features = {
                'test_rna_1': {
                    'conservation': np.random.rand(50),
                    'coevolution': np.random.rand(50, 50),
                    'accessibility': np.random.rand(50),
                    'msa_depth': 10
                },
                'test_rna_2': {
                    'conservation': np.random.rand(30),
                    'coevolution': np.random.rand(30, 30),
                    'accessibility': np.random.rand(30),
                    'msa_depth': 5
                }
            }
            
            # Save features
            output_file = Path(tmpdir) / "test_features.npz"
            np.savez_compressed(output_file, **features)
            
            # Load features
            loaded = np.load(output_file, allow_pickle=True)
            
            assert len(loaded.files) == 2
            assert 'test_rna_1' in loaded.files
            assert 'test_rna_2' in loaded.files
            
            # Check content
            feat1 = loaded['test_rna_1'].item()
            assert feat1['msa_depth'] == 10
            assert feat1['conservation'].shape == (50,)
            
            print(f"  ✓ Saved features for {len(features)} RNAs")
            print(f"  ✓ Loaded features successfully")
            print(f"  ✓ Test 6 passed\n")
            return True
            
        except Exception as e:
            print(f"  ✗ Test 6 failed: {e}\n")
            return False


def main():
    """Run all tests"""
    print("=" * 70)
    print("RNA MSA Pipeline - Test Suite")
    print("=" * 70)
    print()
    
    tests = [
        ("Dependencies", test_dependencies),
        ("RNA Sequence Loading", test_rna_sequence),
        ("Single-Sequence Prediction", test_single_sequence_prediction),
        ("Feature Extraction", test_feature_extraction),
        ("Pipeline Integration", test_pipeline_integration),
        ("Feature Saving/Loading", test_feature_saving),
    ]
    
    results = []
    
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"  ✗ Test failed with exception: {e}\n")
            results.append((test_name, False))
    
    # Summary
    print("=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    print()
    print(f"Total: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n✓ All tests passed! Pipeline is ready to use.")
        print("\nNext steps:")
        print("1. Run pipeline in test mode:")
        print("   python run_rna_msa_pipeline.py --test_mode")
        print("2. Run on full dataset:")
        print("   python run_rna_msa_pipeline.py --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta")
        return 0
    else:
        print(f"\n✗ {total - passed} test(s) failed. Please check the errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
