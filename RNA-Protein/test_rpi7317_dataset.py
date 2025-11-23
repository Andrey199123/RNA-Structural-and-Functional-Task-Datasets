#!/usr/bin/env python3
"""
Test RNA MSA pipeline on RPI7317 dataset
Characterize RNA types, lengths, and MSA generation success
"""

import sys
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from rna_msa.pipeline import RNAMSAPipeline, RNASequence, RNAStrategy

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def analyze_rna_dataset(fasta_path: Path):
    """Analyze RNA dataset characteristics"""
    logger.info(f"Analyzing RNA dataset: {fasta_path}")
    
    rnas = RNASequence.from_fasta(fasta_path)
    
    lengths = [rna.length for rna in rnas]
    
    stats = {
        'total_rnas': len(rnas),
        'min_length': min(lengths),
        'max_length': max(lengths),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'short_rnas': sum(1 for l in lengths if l < 200),
        'medium_rnas': sum(1 for l in lengths if 200 <= l < 1000),
        'long_rnas': sum(1 for l in lengths if l >= 1000),
    }
    
    logger.info("=" * 80)
    logger.info("RNA Dataset Statistics")
    logger.info("=" * 80)
    logger.info(f"Total RNAs: {stats['total_rnas']}")
    logger.info(f"Length range: {stats['min_length']} - {stats['max_length']} nt")
    logger.info(f"Mean length: {stats['mean_length']:.1f} nt")
    logger.info(f"Median length: {stats['median_length']:.1f} nt")
    logger.info(f"")
    logger.info(f"By length category:")
    logger.info(f"  Short (<200nt):     {stats['short_rnas']:4d} ({100*stats['short_rnas']/stats['total_rnas']:.1f}%)")
    logger.info(f"  Medium (200-1000nt): {stats['medium_rnas']:4d} ({100*stats['medium_rnas']/stats['total_rnas']:.1f}%)")
    logger.info(f"  Long (>1000nt):     {stats['long_rnas']:4d} ({100*stats['long_rnas']/stats['total_rnas']:.1f}%)")
    logger.info("=" * 80)
    
    return rnas, stats


def test_msa_pipeline(rnas, pipeline, max_test=None):
    """Test MSA pipeline on RNA dataset"""
    logger.info(f"Testing MSA pipeline on {len(rnas)} RNAs")
    
    if max_test:
        logger.info(f"Limiting to first {max_test} RNAs for testing")
        rnas = rnas[:max_test]
    
    results = {
        'strategy_counts': defaultdict(int),
        'success_counts': defaultdict(int),
        'msa_depths': defaultdict(list),
        'failed_rnas': [],
        'rfam_families': defaultdict(int),
    }
    
    for i, rna in enumerate(rnas, 1):
        logger.info(f"\n[{i}/{len(rnas)}] Processing {rna.id} ({rna.length} nt)")
        
        try:
            result = pipeline.process_rna(rna)
            
            strategy = result.get('strategy', 'unknown')
            success = result.get('success', False)
            msa_depth = result.get('msa_depth', 0)
            
            results['strategy_counts'][strategy] += 1
            
            if success:
                results['success_counts'][strategy] += 1
                results['msa_depths'][strategy].append(msa_depth)
                
                if 'rfam_family' in result:
                    results['rfam_families'][result['rfam_family']] += 1
                
                logger.info(f"  ✓ Success: {strategy}, MSA depth: {msa_depth}")
            else:
                results['failed_rnas'].append((rna.id, rna.length, strategy))
                logger.info(f"  ✗ Failed: {strategy}")
                
        except Exception as e:
            logger.error(f"  ✗ Error processing {rna.id}: {e}")
            results['failed_rnas'].append((rna.id, rna.length, 'error'))
    
    return results


def print_results(results):
    """Print summary of results"""
    logger.info("\n" + "=" * 80)
    logger.info("MSA Pipeline Results")
    logger.info("=" * 80)
    
    # Strategy usage
    logger.info("\nStrategy Usage:")
    for strategy, count in sorted(results['strategy_counts'].items()):
        success = results['success_counts'].get(strategy, 0)
        success_rate = 100 * success / count if count > 0 else 0
        logger.info(f"  {strategy:30s}: {count:4d} RNAs ({success:4d} successful, {success_rate:.1f}%)")
    
    # MSA depths
    logger.info("\nMSA Depth Statistics:")
    for strategy, depths in sorted(results['msa_depths'].items()):
        if depths:
            logger.info(f"  {strategy:30s}: min={min(depths):4d}, max={max(depths):4d}, "
                       f"mean={np.mean(depths):6.1f}, median={np.median(depths):6.1f}")
    
    # Rfam families
    if results['rfam_families']:
        logger.info("\nTop 10 Rfam Families:")
        for family, count in sorted(results['rfam_families'].items(), 
                                    key=lambda x: x[1], reverse=True)[:10]:
            logger.info(f"  {family:20s}: {count:4d} RNAs")
    
    # Failed RNAs
    if results['failed_rnas']:
        logger.info(f"\nFailed RNAs: {len(results['failed_rnas'])}")
        for rna_id, length, strategy in results['failed_rnas'][:10]:
            logger.info(f"  {rna_id:30s} ({length:5d} nt) - {strategy}")
        if len(results['failed_rnas']) > 10:
            logger.info(f"  ... and {len(results['failed_rnas']) - 10} more")
    
    logger.info("=" * 80)


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Test RNA MSA pipeline on RPI7317 dataset')
    parser.add_argument('--rna_fasta', type=Path, 
                       default=Path('data/raw/RPI7317_IncRNA_sequence.fasta'),
                       help='Path to RNA FASTA file')
    parser.add_argument('--rfam_db', type=Path,
                       default=Path('data/databases/Rfam.cm'),
                       help='Path to Rfam.cm database')
    parser.add_argument('--rfam_seed', type=Path,
                       default=Path('data/databases/Rfam.seed'),
                       help='Path to Rfam.seed file')
    parser.add_argument('--output_dir', type=Path,
                       default=Path('data/rna_msa_test'),
                       help='Output directory')
    parser.add_argument('--max_test', type=int, default=None,
                       help='Maximum number of RNAs to test (for quick testing)')
    parser.add_argument('--analyze_only', action='store_true',
                       help='Only analyze dataset, do not run pipeline')
    
    args = parser.parse_args()
    
    # Check input file
    if not args.rna_fasta.exists():
        logger.error(f"RNA FASTA file not found: {args.rna_fasta}")
        return 1
    
    # Analyze dataset
    rnas, stats = analyze_rna_dataset(args.rna_fasta)
    
    if args.analyze_only:
        return 0
    
    # Initialize pipeline
    logger.info("\nInitializing RNA MSA pipeline...")
    
    rfam_db = args.rfam_db if args.rfam_db.exists() else None
    rfam_seed = args.rfam_seed if args.rfam_seed.exists() else None
    
    if rfam_db:
        logger.info(f"  Using Rfam database: {rfam_db}")
    else:
        logger.warning(f"  Rfam database not found: {args.rfam_db}")
    
    if rfam_seed:
        logger.info(f"  Using Rfam.seed: {rfam_seed}")
    else:
        logger.warning(f"  Rfam.seed not found: {args.rfam_seed}")
    
    pipeline = RNAMSAPipeline(
        rfam_db_path=rfam_db,
        rfam_seed_path=rfam_seed,
        output_dir=args.output_dir,
        n_cpu=4
    )
    
    # Test pipeline
    results = test_msa_pipeline(rnas, pipeline, max_test=args.max_test)
    
    # Print results
    print_results(results)
    
    # Save results to JSON
    import json
    results_file = args.output_dir / 'test_results.json'
    
    # Convert to serializable format
    results_json = {
        'dataset_stats': stats,
        'strategy_counts': dict(results['strategy_counts']),
        'success_counts': dict(results['success_counts']),
        'msa_depth_stats': {
            strategy: {
                'min': int(min(depths)),
                'max': int(max(depths)),
                'mean': float(np.mean(depths)),
                'median': float(np.median(depths)),
                'count': len(depths)
            }
            for strategy, depths in results['msa_depths'].items() if depths
        },
        'rfam_families': dict(results['rfam_families']),
        'failed_rnas': [
            {'id': rna_id, 'length': length, 'strategy': strategy}
            for rna_id, length, strategy in results['failed_rnas']
        ]
    }
    
    with open(results_file, 'w') as f:
        json.dump(results_json, f, indent=2)
    
    logger.info(f"\nResults saved to: {results_file}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
