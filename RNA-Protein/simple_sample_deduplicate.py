#!/usr/bin/env python3
"""
Simple approach: Sample 400k sequences, then deduplicate.
Final size will be ~250k-350k after deduplication.
"""

import sys
from pathlib import Path
from Bio import SeqIO
import subprocess
import numpy as np
import argparse
from tqdm import tqdm

CDHIT_IDENTITY = 0.95  # 95% identity - much faster, still removes near-duplicates

def sample_sequences(input_file, output_file, sample_size):
    """Randomly sample sequences from a large FASTA file."""
    print(f"\n{'='*60}")
    print("STEP 1: Random Sampling")
    print(f"{'='*60}")
    print(f"Input: {input_file}")
    print(f"Sample size: {sample_size:,}")
    
    # Count total sequences
    print("\nCounting sequences...")
    total_seqs = sum(1 for _ in SeqIO.parse(input_file, "fasta"))
    print(f"Total sequences: {total_seqs:,}")
    
    if sample_size >= total_seqs:
        print(f"Sample size >= total sequences, using all sequences")
        sample_size = total_seqs
    
    # Generate random indices
    print(f"\nGenerating {sample_size:,} random indices...")
    np.random.seed(42)  # Reproducibility
    selected_indices = set(np.random.choice(total_seqs, sample_size, replace=False))
    
    # Sample sequences
    print(f"Sampling sequences...")
    sampled = []
    for idx, record in enumerate(tqdm(SeqIO.parse(input_file, "fasta"), total=total_seqs)):
        if idx in selected_indices:
            sampled.append(record)
    
    # Write output
    print(f"\nWriting {len(sampled):,} sequences to {output_file}")
    SeqIO.write(sampled, output_file, "fasta")
    
    print(f"\n[DONE] Sampling complete!")
    print(f"  Output: {output_file}")
    print(f"  Sequences: {len(sampled):,}")
    
    return output_file


def deduplicate_sequences(input_file, output_file, num_threads=12):
    """Run CD-HIT-EST deduplication."""
    import re
    import time
    
    print(f"\n{'='*60}")
    print("STEP 2: Deduplication (CD-HIT-EST)")
    print(f"{'='*60}")
    print(f"Input: {input_file}")
    print(f"Identity threshold: {CDHIT_IDENTITY*100:.0f}%")
    print(f"Threads: {num_threads}")
    
    # Find cd-hit-est
    cdhit_path = 'cd-hit-est'
    try:
        result = subprocess.run(['which', 'cd-hit-est'], capture_output=True, text=True)
        if result.returncode == 0:
            cdhit_path = result.stdout.strip()
            print(f"Found cd-hit-est at: {cdhit_path}")
    except:
        pass
    
    # Count input
    input_count = sum(1 for _ in SeqIO.parse(input_file, "fasta"))
    print(f"\nInput sequences: {input_count:,}")
    
    # Run CD-HIT-EST
    cmd = [
        cdhit_path,
        '-i', str(input_file),
        '-o', str(output_file),
        '-c', str(CDHIT_IDENTITY),
        '-n', '10',  # word size for 95% identity
        '-M', '0',
        '-T', str(num_threads),
        '-d', '0',
        '-g', '1',
        '-aS', '0.8',  # alignment coverage for shorter sequence
        '-aL', '0.8'   # alignment coverage for longer sequence
    ]
    
    print(f"\nRunning: {' '.join(cmd)}")
    print("Estimated time: 30 min - 1.5 hours for 400k sequences at 95% identity...")
    print("Progress updates will show below:\n")
    
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        
        start_time = time.time()
        last_progress = 0
        last_time = start_time
        
        for line in process.stdout:
            print(line, end='')
            
            # Track progress and estimate ETA
            match = re.search(r'(\d+)\s+finished', line)
            if match:
                current = int(match.group(1))
                current_time = time.time()
                
                # Calculate progress percentage
                progress = (current / input_count) * 100
                
                # Calculate ETA based on recent progress
                if current > last_progress and current > 1000:
                    elapsed = current_time - start_time
                    rate = current / elapsed  # sequences per second
                    remaining = input_count - current
                    eta_seconds = remaining / rate if rate > 0 else 0
                    
                    # Format ETA
                    eta_hours = int(eta_seconds // 3600)
                    eta_mins = int((eta_seconds % 3600) // 60)
                    
                    print(f"  → Progress: {progress:.1f}% | ETA: {eta_hours}h {eta_mins}m")
                
                last_progress = current
                last_time = current_time
        
        process.wait()
        
        if process.returncode != 0:
            print(f"\n[ERROR] CD-HIT-EST failed with return code {process.returncode}")
            sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] CD-HIT-EST failed: {e}")
        sys.exit(1)
    
    # Count output
    output_count = sum(1 for _ in SeqIO.parse(output_file, "fasta"))
    
    print(f"\n{'='*60}")
    print("[DONE] Deduplication Complete!")
    print(f"{'='*60}")
    print(f"Input sequences:  {input_count:,}")
    print(f"Output sequences: {output_count:,}")
    print(f"Redundancy removed: {(1 - output_count/input_count)*100:.2f}%")
    print(f"Output file: {output_file}")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Simple: Sample 400k, then deduplicate to ~250-350k"
    )
    parser.add_argument('--input', type=str,
                       default='data/rnacentral_subset/02_quality_filtered.fasta',
                       help='Input FASTA file')
    parser.add_argument('--output-dir', type=str,
                       default='data/rnacentral_subset_simple',
                       help='Output directory')
    parser.add_argument('--sample-size', type=int, default=400000,
                       help='Sample size (default: 400k)')
    parser.add_argument('--threads', type=int, default=12,
                       help='Number of CPU threads')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}")
        sys.exit(1)
    
    print(f"\n{'#'*60}")
    print("# Simple RNAcentral Subset Preparation")
    print(f"# Strategy: Sample {args.sample_size:,} → Deduplicate")
    print(f"# Expected final size: ~{int(args.sample_size * 0.7):,}-{int(args.sample_size * 0.85):,}")
    print(f"{'#'*60}")
    
    # Step 1: Sample
    sampled_file = output_dir / "01_sampled.fasta"
    sample_sequences(input_path, sampled_file, args.sample_size)
    
    # Step 2: Deduplicate
    final_file = output_dir / "final_subset.fasta"
    deduplicate_sequences(sampled_file, final_file, args.threads)
    
    print(f"\n{'='*60}")
    print("COMPLETE!")
    print(f"{'='*60}")
    print(f"Final dataset: {final_file}")
    print(f"Next step: Generate statistics and share with Azam")


if __name__ == "__main__":
    main()
