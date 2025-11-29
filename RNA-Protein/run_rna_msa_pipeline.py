#!/usr/bin/env python3
"""
Runs the RNA MSA pipeline on a dataset of RNA sequences.
"""

import argparse
import logging
from pathlib import Path
import sys
import json

sys.path.insert(0, str(Path(__file__).parent / "src"))

from rna_msa.pipeline import RNAMSAPipeline, RNASequence
from rna_msa.msa_utils import check_dependencies, compute_msa_statistics


def setup_logging(verbose: bool = False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('rna_msa_pipeline.log')
        ]
    )


def main():
    parser = argparse.ArgumentParser(
        description="RNA MSA Pipeline for RNA-Protein Interaction Prediction"
    )
    
    parser.add_argument(
        "--rna_fasta",
        type=Path,
        default=Path("RNA-Protein/data/RPI7317/RPI7317_IncRNA_sequence.fasta"),
        help="Path to RNA sequences FASTA file"
    )
    
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("RNA-Protein/data/rna_msa_output"),
        help="Output directory for MSAs and features"
    )
    
    parser.add_argument(
        "--output_features",
        type=Path,
        default=Path("RNA-Protein/data/embeddings/rna_2d_features_rpi7317.npz"),
        help="Output path for 2D feature embeddings"
    )
    
    parser.add_argument(
        "--rfam_db",
        type=Path,
        default=None,
        help="Path to Rfam.cm database (optional)"
    )
    
    parser.add_argument(
        "--rfam_seed",
        type=Path,
        default=None,
        help="Path to Rfam.seed file with family alignments (optional)"
    )
    
    parser.add_argument(
        "--sequence_db",
        type=Path,
        default=None,
        help="Path to sequence database for homology search (optional)"
    )
    
    parser.add_argument(
        "--n_cpu",
        type=int,
        default=4,
        help="Number of CPUs for parallel processing"
    )
    
    parser.add_argument(
        "--check_deps",
        action="store_true",
        help="Check dependencies and exit"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose logging"
    )
    
    parser.add_argument(
        "--test_mode",
        action="store_true",
        help="Test mode: process only first 5 sequences"
    )
    
    parser.add_argument(
        "--allow_fallback",
        action="store_true",
        help="Allow fallback to single-sequence strategy on failure (default: False for research)"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    logger.info("RNA MSA Pipeline for RNA-Protein Interaction Prediction")
    
    # Check dependencies
    logger.info("Checking dependencies...")
    deps = check_dependencies()
    
    for tool, available in deps.items():
        status = "OK" if available else "MISSING"
        logger.info(f"  {status} {tool}")
    
    if args.check_deps:
        return
    
    # Warn about missing tools
    if not deps['infernal']:
        logger.warning("Infernal not found - Rfam search will not work")
    if not deps['RNAfold'] and not deps['python_RNA']:
        logger.warning("ViennaRNA not found - structure prediction will be limited")
    
    # Check input file
    if not args.rna_fasta.exists():
        logger.error(f"RNA FASTA file not found: {args.rna_fasta}")
        return 1
    
    # Initialize pipeline
    logger.info(f"Initializing pipeline...")
    logger.info(f"  RNA FASTA: {args.rna_fasta}")
    logger.info(f"  Output dir: {args.output_dir}")
    logger.info(f"  Rfam DB: {args.rfam_db or 'Not provided'}")
    logger.info(f"  Rfam seed: {args.rfam_seed or 'Not provided'}")
    logger.info(f"  Sequence DB: {args.sequence_db or 'Not provided'}")
    
    pipeline = RNAMSAPipeline(
        rfam_db_path=args.rfam_db,
        rfam_seed_path=args.rfam_seed,
        sequence_db_path=args.sequence_db,
        output_dir=args.output_dir,
        n_cpu=args.n_cpu
    )
    
    # Load RNA sequences
    logger.info(f"Loading RNA sequences from {args.rna_fasta}...")
    sequences = RNASequence.from_fasta(args.rna_fasta)
    logger.info(f"Loaded {len(sequences)} RNA sequences")
    
    if args.test_mode:
        sequences = sequences[:5]
        logger.info(f"TEST MODE: Processing only {len(sequences)} sequences")
    
    lengths = [seq.length for seq in sequences]
    logger.info(f"Sequence length statistics:")
    logger.info(f"  Min: {min(lengths)} nt")
    logger.info(f"  Max: {max(lengths)} nt")
    logger.info(f"  Mean: {sum(lengths)/len(lengths):.1f} nt")
    logger.info(f"  Long RNAs (>1000nt): {sum(1 for l in lengths if l > 1000)}")
    
    logger.info("Processing RNA sequences...")
    
    results = {}
    success_count = 0
    
    for i, rna in enumerate(sequences, 1):
        logger.info(f"\n[{i}/{len(sequences)}] Processing {rna.id} ({rna.length} nt)")
        
        try:
            result = pipeline.process_rna(rna, force_strategy=not args.allow_fallback)
            results[rna.id] = result
            
            if result['success']:
                success_count += 1
                logger.info(f"  Success: {result['strategy']}")
                if 'msa_depth' in result:
                    logger.info(f"    MSA depth: {result['msa_depth']}")
            else:
                logger.warning(f"  Failed: {result.get('error', 'Unknown error')}")
                
        except Exception as e:
            logger.error(f"  Exception: {e}")
            results[rna.id] = {
                'rna_id': rna.id,
                'success': False,
                'error': str(e)
            }
    
    logger.info("Pipeline Summary")
    logger.info(f"Total sequences: {len(sequences)}")
    logger.info(f"Successful: {success_count}")
    logger.info(f"Failed: {len(sequences) - success_count}")
    logger.info(f"Success rate: {success_count/len(sequences)*100:.1f}%")
    
    # Count strategies used
    strategies = {}
    for result in results.values():
        if result['success']:
            strategy = result['strategy']
            strategies[strategy] = strategies.get(strategy, 0) + 1
    
    logger.info("\nStrategies used:")
    for strategy, count in sorted(strategies.items()):
        logger.info(f"  {strategy}: {count}")
    
    summary_file = args.output_dir / "pipeline_results.json"
    with open(summary_file, 'w') as f:
        import numpy as np
        
        def convert_to_json_serializable(obj):
            if isinstance(obj, (int, float, str, bool, type(None))):
                return obj
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_to_json_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_to_json_serializable(item) for item in obj]
            else:
                return str(obj)
        
        json_results = {rna_id: convert_to_json_serializable(result) 
                       for rna_id, result in results.items()}
        
        json.dump(json_results, f, indent=2)
    
    logger.info(f"\nResults summary saved to: {summary_file}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
