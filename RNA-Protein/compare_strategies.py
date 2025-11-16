#!/usr/bin/env python3
"""
Compare MSA Quality Across Different Strategies

This script helps you determine which strategy produces the best MSAs
for your specific dataset.
"""

import sys
import argparse
from pathlib import Path
import json
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "src"))

from rna_msa.pipeline import RNASequence, RNAMSAPipeline, RNAStrategy


def compare_strategies_on_rna(rna, pipeline, strategies):
    """
    Run all available strategies on a single RNA and compare results
    """
    results = {}
    
    for strategy in strategies:
        print(f"\n  Testing {strategy.value}...")
        
        try:
            result = pipeline.process_rna(rna, strategy=strategy, force_strategy=False)
            
            if result['success']:
                msa_depth = result.get('msa_depth', 0)
                
                # Quality assessment
                if msa_depth >= 50:
                    quality = "excellent"
                    score = 1.0
                elif msa_depth >= 20:
                    quality = "good"
                    score = 0.8
                elif msa_depth >= 5:
                    quality = "fair"
                    score = 0.6
                elif msa_depth > 0:
                    quality = "poor"
                    score = 0.4
                else:
                    quality = "single_sequence"
                    score = 0.2
                
                results[strategy.value] = {
                    'success': True,
                    'msa_depth': msa_depth,
                    'quality': quality,
                    'score': score,
                    'strategy_used': result['strategy']
                }
                
                print(f"    ✓ Success: {quality} (depth={msa_depth})")
            else:
                results[strategy.value] = {
                    'success': False,
                    'error': result.get('error', 'unknown')
                }
                print(f"    ✗ Failed: {result.get('error', 'unknown')}")
                
        except Exception as e:
            results[strategy.value] = {
                'success': False,
                'error': str(e)
            }
            print(f"    ✗ Exception: {e}")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Compare MSA generation strategies"
    )
    
    parser.add_argument(
        "--rna_fasta",
        type=Path,
        required=True,
        help="RNA sequences to test"
    )
    
    parser.add_argument(
        "--n_test",
        type=int,
        default=10,
        help="Number of RNAs to test (default: 10)"
    )
    
    parser.add_argument(
        "--rfam_db",
        type=Path,
        default=None,
        help="Path to Rfam database"
    )
    
    parser.add_argument(
        "--sequence_db",
        type=Path,
        default=None,
        help="Path to sequence database"
    )
    
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("data/strategy_comparison"),
        help="Output directory"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("RNA MSA Strategy Comparison")
    print("=" * 70)
    print()
    
    # Check what's available
    print("Available resources:")
    print(f"  Rfam DB: {'✓' if args.rfam_db and args.rfam_db.exists() else '✗'}")
    print(f"  Sequence DB: {'✓' if args.sequence_db and args.sequence_db.exists() else '✗'}")
    print()
    
    # Load RNAs
    print(f"Loading RNAs from {args.rna_fasta}...")
    sequences = RNASequence.from_fasta(args.rna_fasta)
    print(f"Loaded {len(sequences)} RNAs")
    
    # Select test set
    test_sequences = sequences[:args.n_test]
    print(f"Testing on {len(test_sequences)} RNAs")
    print()
    
    # Initialize pipeline
    pipeline = RNAMSAPipeline(
        rfam_db_path=args.rfam_db,
        sequence_db_path=args.sequence_db,
        output_dir=args.output_dir,
        n_cpu=4
    )
    
    # Determine which strategies to test
    strategies = [RNAStrategy.SINGLE_SEQUENCE]  # Always available
    
    if args.rfam_db and args.rfam_db.exists():
        strategies.append(RNAStrategy.INFERNAL_RFAM)
    
    if args.sequence_db and args.sequence_db.exists():
        strategies.append(RNAStrategy.RNACMAP_STYLE)
    
    print(f"Testing {len(strategies)} strategies:")
    for s in strategies:
        print(f"  - {s.value}")
    print()
    
    # Test each RNA
    all_results = {}
    
    for i, rna in enumerate(test_sequences, 1):
        print(f"[{i}/{len(test_sequences)}] {rna.id} ({rna.length} nt)")
        
        results = compare_strategies_on_rna(rna, pipeline, strategies)
        all_results[rna.id] = results
    
    # Analyze results
    print("\n" + "=" * 70)
    print("Strategy Comparison Summary")
    print("=" * 70)
    print()
    
    strategy_stats = {s.value: {'success': 0, 'total_depth': 0, 'scores': []} 
                      for s in strategies}
    
    for rna_id, results in all_results.items():
        for strategy, result in results.items():
            if result['success']:
                strategy_stats[strategy]['success'] += 1
                strategy_stats[strategy]['total_depth'] += result.get('msa_depth', 0)
                strategy_stats[strategy]['scores'].append(result.get('score', 0))
    
    # Print summary
    for strategy in strategies:
        stats = strategy_stats[strategy.value]
        n_success = stats['success']
        n_total = len(test_sequences)
        
        print(f"\n{strategy.value}:")
        print(f"  Success rate: {n_success}/{n_total} ({n_success/n_total*100:.1f}%)")
        
        if n_success > 0:
            avg_depth = stats['total_depth'] / n_success
            avg_score = np.mean(stats['scores'])
            print(f"  Average MSA depth: {avg_depth:.1f}")
            print(f"  Average quality score: {avg_score:.2f}")
    
    # Save results
    output_file = args.output_dir / "strategy_comparison.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\nDetailed results saved to: {output_file}")
    
    # Recommendation
    print("\n" + "=" * 70)
    print("Recommendation")
    print("=" * 70)
    
    best_strategy = max(strategies, 
                       key=lambda s: np.mean(strategy_stats[s.value]['scores']) 
                       if strategy_stats[s.value]['scores'] else 0)
    
    print(f"\nBest strategy for your dataset: {best_strategy.value}")
    print()
    print("However, the pipeline automatically selects the best strategy")
    print("for each RNA, so you don't need to choose manually.")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
