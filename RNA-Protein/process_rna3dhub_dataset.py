#!/usr/bin/env python3
"""
Process RNA 3D Hub dataset (13,782 RNAs) using hybrid MSA generation strategy.
Generates A3M files compatible with trRosettaRNA2 feature extraction.

Fixes implemented:
1. BLAST+MAFFT: Query sequence included in MAFFT input, alignment columns preserved
2. Infernal+Rfam: Stockholm alignment columns preserved during conversion
3. Ambiguity codes (WRYKMSBDHVN) normalized to 'N'
4. Each sequence written on ONE line
5. Query sequence is first in A3M file with NO gaps
6. All sequences have same length (alignment preserved)
"""

import subprocess
from pathlib import Path
import time
import sys
import csv
import argparse
from typing import Dict, List

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))
from rna_msa.alignment_formats import fasta_to_a3m, stockholm_to_a3m, validate_a3m


def run_infernal_rfam(rna: Dict, output_dir: Path, rfam_cm: Path, rfam_seed: Path, 
                      cmscan_path: str, timeout: int = 60) -> tuple:
    """
    Try Infernal+Rfam approach.
    
    Returns:
        (success: bool, a3m_file: Path or None, n_sequences: int, error: str or None)
    """
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    try:
        # Run cmscan
        infernal_out = output_dir / f"{rna['id']}_infernal.tblout"
        cmd = [
            cmscan_path,
            '--tblout', str(infernal_out),
            '--cut_ga',
            '--nohmmonly',
            '--cpu', '4',
            str(rfam_cm),
            str(query_file)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        
        # Parse results
        families_found = []
        with open(infernal_out, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 16:
                    family_acc = parts[1]
                    families_found.append(family_acc)
        
        if not families_found:
            return False, None, 0, "No Rfam family found"
        
        best_acc = families_found[0]
        
        # Extract seed alignment from Rfam.seed
        sto_file = output_dir / f"{rna['id']}_{best_acc}.sto"
        with open(rfam_seed, 'r', errors='ignore') as fin:
            with open(sto_file, 'w') as fout:
                in_family = False
                for line in fin:
                    if line.startswith(f'#=GF AC   {best_acc}'):
                        in_family = True
                        fout.write('# STOCKHOLM 1.0\n')
                    if in_family:
                        fout.write(line)
                        if line.startswith('//'):
                            break
        
        # CRITICAL FIX: Convert Stockholm to A3M preserving alignment columns
        # Query sequence is provided separately since it's not in seed alignment
        a3m_file = output_dir / f"{rna['id']}_{best_acc}.a3m"
        stockholm_to_a3m(sto_file, a3m_file, query_id=rna['id'], 
                        query_sequence=rna['seq'])
        
        # Validate
        stats = validate_a3m(a3m_file)
        if stats['valid']:
            return True, a3m_file, stats['n_sequences'], None
        else:
            return False, None, 0, f"Invalid A3M: {stats.get('error', 'Unknown')}"
    
    except subprocess.TimeoutExpired:
        return False, None, 0, "Timeout"
    except Exception as e:
        return False, None, 0, str(e)


def run_blast_mafft(rna: Dict, output_dir: Path, rnacentral_db: Path,
                    blastn_path: str, blastdbcmd_path: str, mafft_path: str,
                    min_identity: float = 70.0, min_homologs: int = 10,
                    timeout: int = 180) -> tuple:
    """
    Try BLAST+MAFFT approach.
    
    Returns:
        (success: bool, a3m_file: Path or None, n_sequences: int, error: str or None)
    """
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    try:
        # Run BLAST
        blast_out = output_dir / f"{rna['id']}_blast.txt"
        cmd = [
            blastn_path,
            '-query', str(query_file),
            '-db', str(rnacentral_db),
            '-evalue', '0.001',
            '-outfmt', '6 sseqid pident length',
            '-max_target_seqs', '100',
            '-num_threads', '4'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        # Parse homologs
        homolog_ids = set()
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 3:
                    seq_id = parts[0]
                    pident = float(parts[1])
                    if pident >= min_identity:
                        homolog_ids.add(seq_id)
        
        n_homologs = len(homolog_ids)
        
        if n_homologs < min_homologs:
            return False, None, 0, f"Insufficient homologs ({n_homologs} < {min_homologs})"
        
        # Extract homologs
        id_file = output_dir / f"{rna['id']}_ids.txt"
        with open(id_file, 'w') as f:
            for seq_id in homolog_ids:
                f.write(f"{seq_id}\n")
        
        homologs_file = output_dir / f"{rna['id']}_homologs.fa"
        cmd = [
            blastdbcmd_path,
            '-db', str(rnacentral_db),
            '-entry_batch', str(id_file),
            '-line_length', '100000',
            '-out', str(homologs_file)
        ]
        subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        # CRITICAL FIX: Add query sequence to homologs file for MAFFT
        # This ensures alignment is anchored to the query
        with open(homologs_file, 'r') as f:
            homolog_content = f.read()
        
        with open(homologs_file, 'w') as f:
            # Write query first
            f.write(f">{rna['id']}\n{rna['seq']}\n")
            # Then write homologs
            f.write(homolog_content)
        
        # Run MAFFT - output is aligned FASTA with gaps
        aligned_file = output_dir / f"{rna['id']}_aligned.fa"
        cmd = [mafft_path, '--auto', '--thread', '4', str(homologs_file)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        
        with open(aligned_file, 'w') as f:
            f.write(result.stdout)
        
        # CRITICAL FIX: Convert aligned FASTA to A3M preserving alignment columns
        # Do NOT remove gaps - they are essential for trRosettaRNA
        a3m_file = output_dir / f"{rna['id']}_blast.a3m"
        fasta_to_a3m(aligned_file, a3m_file, query_id=rna['id'])
        
        # Validate
        stats = validate_a3m(a3m_file)
        if stats['valid']:
            return True, a3m_file, stats['n_sequences'], None
        else:
            return False, None, 0, f"Invalid A3M: {stats.get('error', 'Unknown')}"
    
    except subprocess.TimeoutExpired:
        return False, None, 0, "Timeout"
    except Exception as e:
        return False, None, 0, str(e)


def process_rna3dhub_dataset(
    input_fasta: Path,
    output_dir: Path,
    rfam_cm: Path,
    rfam_seed: Path,
    rnacentral_db: Path,
    cmscan_path: str,
    blastn_path: str,
    blastdbcmd_path: str,
    mafft_path: str,
    start_idx: int = 0,
    end_idx: int = None
):
    """
    Process RNA 3D Hub dataset using hybrid MSA generation strategy.
    
    Args:
        input_fasta: FASTA file with RNA sequences from RNA 3D Hub
        output_dir: Directory for output A3M files
        rfam_cm: Path to Rfam.cm covariance model file
        rfam_seed: Path to Rfam.seed alignment file
        rnacentral_db: Path to RNAcentral BLAST database
        cmscan_path: Path to cmscan executable
        blastn_path: Path to blastn executable
        blastdbcmd_path: Path to blastdbcmd executable
        mafft_path: Path to mafft executable
        start_idx: Start index for processing (for parallelization)
        end_idx: End index for processing (for parallelization)
    """
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Parse input FASTA
    print("Loading RNA sequences...")
    rnas = []
    current_id = None
    current_seq = []
    
    with open(input_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id and current_seq:
                    seq = ''.join(current_seq)
                    rnas.append({
                        'id': current_id,
                        'seq': seq,
                        'length': len(seq)
                    })
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id and current_seq:
            seq = ''.join(current_seq)
            rnas.append({
                'id': current_id,
                'seq': seq,
                'length': len(seq)
            })
    
    print(f"Loaded {len(rnas)} RNA sequences")
    
    # Apply start/end indices
    if end_idx is None:
        end_idx = len(rnas)
    rnas = rnas[start_idx:end_idx]
    print(f"Processing RNAs {start_idx} to {end_idx} ({len(rnas)} total)")
    
    # Process each RNA
    results = []
    infernal_success = 0
    blast_success = 0
    total_success = 0
    
    for i, rna in enumerate(rnas, start_idx + 1):
        print(f"\n[{i}/{end_idx}] {rna['id']} ({rna['length']} nt)")
        
        start_time = time.time()
        
        # STEP 1: Try Infernal+Rfam
        print(f"  [1] Trying Infernal+Rfam...", end='', flush=True)
        success, a3m_file, n_seqs, error = run_infernal_rfam(
            rna, output_dir, rfam_cm, rfam_seed, cmscan_path
        )
        
        if success:
            infernal_success += 1
            total_success += 1
            final_method = 'Infernal'
            print(f" ✓ Success ({n_seqs} seqs)")
        else:
            print(f" ✗ {error}")
            
            # STEP 2: Try BLAST+MAFFT
            print(f"  [2] Trying BLAST+MAFFT...", end='', flush=True)
            success, a3m_file, n_seqs, error = run_blast_mafft(
                rna, output_dir, rnacentral_db, blastn_path, 
                blastdbcmd_path, mafft_path
            )
            
            if success:
                blast_success += 1
                total_success += 1
                final_method = 'BLAST'
                print(f" ✓ Success ({n_seqs} seqs)")
            else:
                print(f" ✗ {error}")
                final_method = None
        
        elapsed = time.time() - start_time
        
        if not success:
            print(f"  [FINAL] ✗ Both methods failed")
        
        print(f"  Time: {elapsed:.1f}s")
        
        results.append({
            'rna_id': rna['id'],
            'length': rna['length'],
            'infernal_success': final_method == 'Infernal',
            'blast_success': final_method == 'BLAST',
            'hybrid_success': success,
            'final_method': final_method,
            'a3m_file': str(a3m_file.name) if a3m_file else None,
            'n_sequences': n_seqs if success else 0,
            'time': elapsed
        })
        
        # Save intermediate results every 100 RNAs
        if i % 100 == 0:
            tsv_file = output_dir / f'rna3dhub_results_{start_idx}_{end_idx}.tsv'
            with open(tsv_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=[
                    'rna_id', 'length', 'infernal_success', 'blast_success', 
                    'hybrid_success', 'final_method', 'a3m_file', 'n_sequences', 'time'
                ], delimiter='\t')
                writer.writeheader()
                writer.writerows(results)
    
    # Final summary
    print(f"\n{'='*60}")
    print(f"SUMMARY (RNAs {start_idx}-{end_idx})")
    print(f"{'='*60}")
    print(f"Total RNAs processed: {len(rnas)}")
    print(f"Infernal success: {infernal_success}/{len(rnas)} ({infernal_success/len(rnas)*100:.1f}%)")
    print(f"BLAST success: {blast_success}/{len(rnas)} ({blast_success/len(rnas)*100:.1f}%)")
    print(f"HYBRID total success: {total_success}/{len(rnas)} ({total_success/len(rnas)*100:.1f}%)")
    print(f"Failed: {len(rnas) - total_success}/{len(rnas)} ({(len(rnas) - total_success)/len(rnas)*100:.1f}%)")
    
    # Save final results
    tsv_file = output_dir / f'rna3dhub_results_{start_idx}_{end_idx}.tsv'
    with open(tsv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'rna_id', 'length', 'infernal_success', 'blast_success', 
            'hybrid_success', 'final_method', 'a3m_file', 'n_sequences', 'time'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults saved to: {tsv_file}")
    print(f"A3M files saved to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process RNA 3D Hub dataset with hybrid MSA generation"
    )
    parser.add_argument(
        '--input', '-i', type=Path, required=True,
        help='Input FASTA file with RNA sequences'
    )
    parser.add_argument(
        '--output', '-o', type=Path, required=True,
        help='Output directory for A3M files'
    )
    parser.add_argument(
        '--rfam-cm', type=Path, required=True,
        help='Path to Rfam.cm covariance model file'
    )
    parser.add_argument(
        '--rfam-seed', type=Path, required=True,
        help='Path to Rfam.seed alignment file'
    )
    parser.add_argument(
        '--rnacentral-db', type=Path, required=True,
        help='Path to RNAcentral BLAST database'
    )
    parser.add_argument(
        '--cmscan', type=str, default='cmscan',
        help='Path to cmscan executable'
    )
    parser.add_argument(
        '--blastn', type=str, default='blastn',
        help='Path to blastn executable'
    )
    parser.add_argument(
        '--blastdbcmd', type=str, default='blastdbcmd',
        help='Path to blastdbcmd executable'
    )
    parser.add_argument(
        '--mafft', type=str, default='mafft',
        help='Path to mafft executable'
    )
    parser.add_argument(
        '--start', type=int, default=0,
        help='Start index for processing (for parallelization)'
    )
    parser.add_argument(
        '--end', type=int, default=None,
        help='End index for processing (for parallelization)'
    )
    
    args = parser.parse_args()
    
    process_rna3dhub_dataset(
        input_fasta=args.input,
        output_dir=args.output,
        rfam_cm=args.rfam_cm,
        rfam_seed=args.rfam_seed,
        rnacentral_db=args.rnacentral_db,
        cmscan_path=args.cmscan,
        blastn_path=args.blastn,
        blastdbcmd_path=args.blastdbcmd,
        mafft_path=args.mafft,
        start_idx=args.start,
        end_idx=args.end
    )
