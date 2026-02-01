#!/usr/bin/env python3
"""
RNA structure prediction - Phase 2 only

Skips MSA generation, directly predicts structures from existing .a3m files.
Includes length filtering to skip sequences that are too long for trRosettaRNA2.

Usage:
    conda activate rna_msa
    python generate_rna_structures_phase2.py --resume
    python generate_rna_structures_phase2.py --resume --max-length 1500
"""

import argparse
import subprocess
import json
import os
from pathlib import Path
from datetime import datetime
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp


TRROSETTA_DIR = Path('trRosettaRNA2').resolve()


def get_sequence_info(a3m_file, rna_id):
    """Get query sequence length and validate MSA has query sequence and is aligned"""
    with open(a3m_file) as f:
        lines = f.readlines()
    
    if not lines:
        return None, "empty_file"
    
    # Check if first sequence is the query
    first_header = lines[0].strip() if lines else ""
    if not first_header.startswith(f'>{rna_id}'):
        return None, "missing_query"
    
    # Parse all sequences and check alignment
    sequences = []
    current_seq = ""
    for line in lines:
        if line.startswith('>'):
            if current_seq:
                sequences.append(current_seq)
            current_seq = ""
        else:
            current_seq += line.strip()
    if current_seq:
        sequences.append(current_seq)
    
    if not sequences:
        return None, "empty_sequence"
    
    # Check if all sequences have the same length (aligned)
    seq_lengths = set(len(s) for s in sequences)
    if len(seq_lengths) > 1:
        return None, "unaligned"
    
    query_len = len(sequences[0].replace('-', '').replace('.', ''))
    return query_len, "ok"


def worker_init(gpu_queue):
    """Initialize worker with a GPU from the queue"""
    global WORKER_GPU
    WORKER_GPU = gpu_queue.get()


def predict_structure(args_tuple):
    """Predict structure for a single RNA"""
    rna_id, a3m_file, pdb_dir, timeout = args_tuple
    pdb_dir = Path(pdb_dir)
    pdb_file = pdb_dir / f"{rna_id}.pdb"
    
    # Get GPU assigned to this worker
    gpu_id = globals().get('WORKER_GPU', None)
    
    if pdb_file.exists():
        return rna_id, True, "cached"
    
    try:
        # Resolve paths relative to current working directory
        trrosetta_dir = Path('trRosettaRNA2').resolve()
        abs_a3m = Path(a3m_file).resolve()
        abs_out = (pdb_dir / rna_id).resolve()
        
        # Use relative path for params when running from trRosettaRNA2 dir
        cmd = [
            'python', '-m', 'trRNA2.predict',
            '-mdir', 'params',  # Relative to cwd (trRosettaRNA2)
            '-i', str(abs_a3m),
            '-o', str(abs_out),
            '-relax_steps', '0',
            '-cpu', '2'  # Reduced for parallel runs
        ]
        
        # Set environment to use specific GPU (or CPU if gpu_id is None)
        env = dict(os.environ)
        if gpu_id is not None:
            env['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
        else:
            env['CUDA_VISIBLE_DEVICES'] = ''
        
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True,
            timeout=timeout, 
            cwd=str(trrosetta_dir),
            env=env
        )
        
        predicted_pdb = pdb_dir / rna_id / "model_1_unrelaxed.pdb"
        if predicted_pdb.exists():
            pdb_file.write_text(predicted_pdb.read_text())
            return rna_id, True, "predicted"
        else:
            # Check for error in output
            if "out of memory" in result.stderr.lower():
                return rna_id, False, "OOM"
            return rna_id, False, "no_output"
            
    except subprocess.TimeoutExpired:
        return rna_id, False, "timeout"
    except Exception as e:
        return rna_id, False, str(e)[:50]


def main():
    parser = argparse.ArgumentParser(description="Phase 2: Structure prediction from existing MSAs")
    parser.add_argument('--resume', action='store_true', help="Resume from progress file")
    parser.add_argument('--max-length', type=int, default=2000,
                       help="Max sequence length to attempt (default: 2000)")
    parser.add_argument('--timeout', type=int, default=1800,
                       help="Timeout per prediction in seconds (default: 1800)")
    parser.add_argument('--workers', type=int, default=1,
                       help="Number of parallel workers (default: 1)")
    parser.add_argument('--gpus', type=str, default=None,
                       help="Comma-separated GPU IDs to use (e.g., '0,1,2,3'). If not set, uses CPU.")
    args = parser.parse_args()
    
    # Parse GPU list
    gpu_list = None
    if args.gpus:
        gpu_list = [int(g.strip()) for g in args.gpus.split(',')]
    
    print("=" * 80)
    print("RNA Structure Prediction - Phase 2 Only")
    print("=" * 80)
    print(f"Max sequence length: {args.max_length} nt")
    print(f"Timeout per RNA: {args.timeout}s")
    print(f"Parallel workers: {args.workers}")
    if gpu_list:
        print(f"Using GPUs: {gpu_list}")
    else:
        print("Using CPU only")
    print()
    
    # Paths
    output_dir = Path("RNA-Protein/data/structural_features")
    msa_dir = output_dir / "msas"
    pdb_dir = output_dir / "rna_structures"
    progress_file = output_dir / "rna_generation_progress.json"
    
    pdb_dir.mkdir(parents=True, exist_ok=True)
    
    # Load progress
    completed_structures = set()
    if args.resume and progress_file.exists():
        with open(progress_file) as f:
            progress = json.load(f)
            completed_structures = set(progress.get("completed", []))
        print(f"Already completed: {len(completed_structures)} structures")
    
    # Find all MSA files
    a3m_files = list(msa_dir.glob("*.a3m"))
    print(f"Found {len(a3m_files)} MSA files")
    
    # Filter and categorize
    to_process = []
    skipped_completed = 0
    skipped_too_long = 0
    skipped_bad_msa = 0
    
    for a3m_file in a3m_files:
        rna_id = a3m_file.stem
        
        if rna_id in completed_structures:
            skipped_completed += 1
            continue
        
        # Check if already has PDB
        if (pdb_dir / f"{rna_id}.pdb").exists():
            completed_structures.add(rna_id)
            skipped_completed += 1
            continue
        
        seq_len, status = get_sequence_info(a3m_file, rna_id)
        
        if status == "unaligned":
            skipped_bad_msa += 1
            continue
        elif status != "ok":
            skipped_bad_msa += 1
            continue
            
        if seq_len > args.max_length:
            skipped_too_long += 1
            continue
        
        to_process.append((rna_id, str(a3m_file), seq_len))
    
    # Sort by length (shorter first - more likely to succeed)
    to_process.sort(key=lambda x: x[2])
    
    print(f"Skipped (already done): {skipped_completed}")
    print(f"Skipped (bad MSA - missing query): {skipped_bad_msa}")
    print(f"Skipped (too long >{args.max_length}nt): {skipped_too_long}")
    print(f"To process: {len(to_process)}")
    print()
    
    if not to_process:
        print("Nothing to process!")
        return
    
    # Process
    print("=" * 80)
    print("Structure Prediction")
    print("=" * 80)
    
    success_count = 0
    fail_count = 0
    start_time = time.time()
    processed = 0
    
    if args.workers == 1:
        # Sequential processing
        # Set global GPU for single worker
        global WORKER_GPU
        WORKER_GPU = gpu_list[0] if gpu_list else None
        
        for i, (rna_id, a3m_file, seq_len) in enumerate(to_process, 1):
            elapsed = time.time() - start_time
            rate = i / elapsed if elapsed > 0 else 0
            remaining = len(to_process) - i
            eta = remaining / rate / 3600 if rate > 0 else 0
            
            print(f"[{i}/{len(to_process)}] {rna_id} ({seq_len}nt)...", end='', flush=True)
            
            _, success, reason = predict_structure((rna_id, a3m_file, str(pdb_dir), args.timeout))
            
            if success:
                success_count += 1
                completed_structures.add(rna_id)
                print(f" ✓ {reason} | {success_count}/{i} | ETA: {eta:.1f}h")
                
                # Save progress
                with open(progress_file, 'w') as f:
                    json.dump({
                        "completed": list(completed_structures),
                        "last_updated": datetime.now().isoformat()
                    }, f)
            else:
                fail_count += 1
                print(f" ✗ {reason}")
    else:
        # Parallel processing - each worker gets a dedicated GPU
        work_items = [(rna_id, a3m_file, str(pdb_dir), args.timeout) 
                      for rna_id, a3m_file, seq_len in to_process]
        
        # Create a queue with GPU IDs for workers
        gpu_queue = mp.Manager().Queue()
        if gpu_list:
            # Assign GPUs to workers (cycle if more workers than GPUs)
            for i in range(args.workers):
                gpu_queue.put(gpu_list[i % len(gpu_list)])
        else:
            # CPU only - put None for each worker
            for i in range(args.workers):
                gpu_queue.put(None)
        
        with ProcessPoolExecutor(max_workers=args.workers, 
                                  initializer=worker_init, 
                                  initargs=(gpu_queue,)) as executor:
            futures = {executor.submit(predict_structure, item): item for item in work_items}
            
            for future in as_completed(futures):
                processed += 1
                rna_id, success, reason = future.result()
                
                # Find seq_len for this rna_id
                seq_len = next(s for r, _, s in to_process if r == rna_id)
                
                elapsed = time.time() - start_time
                rate = processed / elapsed if elapsed > 0 else 0
                remaining = len(to_process) - processed
                eta = remaining / rate / 3600 if rate > 0 else 0
                
                if success:
                    success_count += 1
                    completed_structures.add(rna_id)
                    print(f"[{processed}/{len(to_process)}] {rna_id} ({seq_len}nt)... ✓ {reason} | {success_count}/{processed} | ETA: {eta:.1f}h")
                    
                    # Save progress
                    with open(progress_file, 'w') as f:
                        json.dump({
                            "completed": list(completed_structures),
                            "last_updated": datetime.now().isoformat()
                        }, f)
                else:
                    fail_count += 1
                    print(f"[{processed}/{len(to_process)}] {rna_id} ({seq_len}nt)... ✗ {reason}")
    
    # Summary
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"This run: {success_count} succeeded, {fail_count} failed")
    print(f"Total completed: {len(completed_structures)}")
    print()


if __name__ == "__main__":
    main()
