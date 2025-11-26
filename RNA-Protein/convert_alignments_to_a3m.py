#!/usr/bin/env python3
"""
Convert all Stockholm alignments to A3M format
Run this after MSA generation to get .a3m files
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from rna_msa.alignment_formats import stockholm_to_a3m, validate_a3m
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def convert_directory(input_dir: Path, output_dir: Path = None):
    """
    Convert all .sto files in directory to .a3m format
    
    Args:
        input_dir: Directory containing .sto files
        output_dir: Output directory (if None, uses same as input)
    """
    if output_dir is None:
        output_dir = input_dir
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all Stockholm files
    sto_files = list(input_dir.glob("*.sto")) + list(input_dir.glob("**/*.sto"))
    
    if not sto_files:
        logger.warning(f"No .sto files found in {input_dir}")
        return
    
    logger.info(f"Found {len(sto_files)} Stockholm files")
    
    converted = 0
    failed = 0
    
    for sto_file in sto_files:
        # Create output filename
        a3m_file = output_dir / sto_file.with_suffix('.a3m').name
        
        try:
            logger.info(f"Converting {sto_file.name}...")
            stockholm_to_a3m(sto_file, a3m_file)
            
            # Validate
            stats = validate_a3m(a3m_file)
            if stats['valid']:
                logger.info(f"  ✓ {stats['n_sequences']} sequences, query length {stats['query_length']}")
                converted += 1
            else:
                logger.error(f"  ✗ Invalid A3M: {stats.get('error', 'Unknown error')}")
                failed += 1
                
        except Exception as e:
            logger.error(f"  ✗ Failed: {e}")
            failed += 1
    
    logger.info(f"\nConversion complete:")
    logger.info(f"  Converted: {converted}")
    logger.info(f"  Failed: {failed}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Convert Stockholm alignments to A3M format")
    parser.add_argument("input_dir", type=Path, help="Directory containing .sto files")
    parser.add_argument("--output_dir", type=Path, help="Output directory (default: same as input)")
    
    args = parser.parse_args()
    
    if not args.input_dir.exists():
        logger.error(f"Input directory not found: {args.input_dir}")
        sys.exit(1)
    
    convert_directory(args.input_dir, args.output_dir)
