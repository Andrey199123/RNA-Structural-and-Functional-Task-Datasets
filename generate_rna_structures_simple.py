#!/usr/bin/env python3
"""
Simple RNA structure generation for RPI7317 - uses direct subprocess calls

Usage:
    conda activate rna_msa
    python generate_rna_structures_simple.py --workers 12
"""

import argparse
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO
import json
from datetime import datetime
import time
import multiprocessing as mp

# Paths
RFAM_CM = Path('RNA-Protein/data/rfam/Rfam.cm')
RFAM_SEED = Path('RNA-Protein/data/rfam/Rfam.seed')
RNACENTRAL_DB = Path('RNA-Protein/data/rnacentral/rnacentral_active')
TRROSETTA_DIR = Path('trRosettaRNA2')

# Add src to path for alignment format conversion
sys.path.insert(0, 'RNA-Protein/src')
from rna_msa.alignment_formats import stockholm_to_a3m, fasta_to_a3m, validate_a3m


def process_rna(args_tuple):
    """Process a single RNA"""
    rna_id, sequence, output_dir = args_tuple
    
    msa_dir = output_dir / "msas"
    pdb_dir = output_dir / "rna_structures"
    
    result = {
        "rna_id": rna_id,
        "msa_success": False,
        "structure_success": False,
        "method": None
    }
    
    # Check if already done
    pdb_file = pdb_dir / f"{rna_id}.pdb"
    if pdb_file.exists():
        result["msa_success"] = True
        result["structure_success"] = True
        result["method"] = "cached"
        return result
    
    # Write query
    query_file = msa_dir / f"{rna_id}_query.fa"
    query_file.write_text(f">{rna_id}\n{sequence}\n")
    
    a3m_file = None
    
    # Try Infernal first
    try:
        infernal_out = msa_dir / f"{rna_id}_infernal.tblout"
        cmd = [
            'cmscan',
            '--tblout', str(infernal_out),
            '--cut_ga',
            '--nohmmonly',
            '--cpu', '1',
            str(RFAM_CM),
            str(query_file)
        ]
        
        subprocess.run(cmd, capture_output=True, timeout=60)
        
        # Parse for families
        families = []
        with open(infernal_out) as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 16:
                        families.append(parts[1])
        
        if families:
            best_acc = families[0]
            
            # Extract seed alignment
            sto_file = msa_dir / f"{rna_id}_{best_acc}.sto"
            with open(RFAM_SEED, errors='ignore') as fin:
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
            
            # Convert to A3M
            a3m_file = msa_dir / f"{rna_id}.a3m"
            stockholm_to_a3m(sto_file, a3m_file, query_id=None)
            
            stats = validate_a3m(a3m_file)
            if stats['valid']:
                result["msa_success"] = True
                result["method"] = "infernal"
    except Exception:
        pass
    
    # Try BLAST if Infernal failed
    if not result["msa_success"]:
        try:
            blast_out = msa_dir / f"{rna_id}_blast.txt"
            cmd = [
                'blastn',
                '-query', str(query_file),
                '-db', str(RNACENTRAL_DB),
                '-evalue', '0.001',
                '-outfmt', '6 sseqid pident',
                '-max_target_seqs', '100',
                '-num_threads', '1'
            ]
            
            blast_result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            # Get homologs
            homolog_ids = set()
            for line in blast_result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 2 and float(parts[1]) >= 70:
                        homolog_ids.add(parts[0])
            
            if len(homolog_ids) >= 10:
                # Extract sequences
                id_file = msa_dir / f"{rna_id}_ids.txt"
                id_file.write_text('\n'.join(homolog_ids))
                
                homologs_file = msa_dir / f"{rna_id}_homologs.fa"
                cmd = [
                    'blastdbcmd',
                    '-db', str(RNACENTRAL_DB),
                    '-entry_batch', str(id_file),
                    '-line_length', '100000',
                    '-out', str(homologs_file)
                ]
                subprocess.run(cmd, capture_output=True, timeout=60)
                
                # Add query sequence to homologs
                with open(homologs_file, 'r') as f:
                    homolog_content = f.read()
                with open(homologs_file, 'w') as f:
                    f.write(f">{rna_id}\n{sequence}\n")
                    f.write(homolog_content)
                
                # Align with MAFFT
                cmd = ['mafft', '--auto', '--thread', '1', str(homologs_file)]
                mafft_result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                
                aligned_file = msa_dir / f"{rna_id}_aligned.fa"
                aligned_file.write_text(mafft_result.stdout)
                
                # Convert to A3M
                a3m_file = msa_dir / f"{rna_id}.a3m"
                fasta_to_a3m(aligned_file, a3m_file, query_id=rna_id)
                
                stats = validate_a3m(a3m_file)
                if stats['valid']:
                    result["msa_success"] = True
                    result["method"] = "blast"
        except Exception:
            pass
    
    # Predict structure if MSA succeeded
    if result["msa_success"] and a3m_file and a3m_file.exists():
        try:
            trrosetta_dir = Path("trRosettaRNA2").resolve()
            abs_a3m = a3m_file.resolve()
            abs_out = (pdb_dir / rna_id).resolve()
            abs_params = (trrosetta_dir / 'params').resolve()
            
            cmd = [
                'python', '-m', 'trRNA2.predict',
                '-mdir', str(abs_params),
                '-i', str(abs_a3m),
                '-o', str(abs_out),
                '-relax_steps', '0'
            ]
            
            subprocess.run(cmd, capture_output=True, timeout=1800, cwd=str(trrosetta_dir))
            
            # trRosettaRNA2 outputs to model_1_unrelaxed.pdb
            predicted_pdb = pdb_dir / rna_id / "model_1_unrelaxed.pdb"
            if predicted_pdb.exists():
                pdb_file.write_text(predicted_pdb.read_text())
                result["structure_success"] = True
        except Exception:
            pass
    
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=mp.cpu_count())
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    
    print("=" * 80)
    print("RPI7317 RNA Structure Generation (Simple Version)")
    print("=" * 80)
    print()
    
    # Setup directories
    output_dir = Path("RNA-Protein/data/structural_features")
    msa_dir = output_dir / "msas"
    pdb_dir = output_dir / "rna_structures"
    progress_file = output_dir / "rna_generation_progress.json"
    
    msa_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)
    
    # Load progress
    completed = set()
    if args.resume and progress_file.exists():
        with open(progress_file) as f:
            progress = json.load(f)
            completed = set(progress.get("completed", []))
        print(f"RESUME: Skipping {len(completed)} completed RNAs")
        print()
    
    # Load RNAs
    print("Loading RNAs...")
    rnas = {}
    for record in SeqIO.parse("RNA-Protein/data/RPI7317/RPI7317_IncRNA_sequence.fasta", "fasta"):
        if record.id not in completed:
            rnas[record.id] = str(record.seq)
    
    print(f"Processing {len(rnas)} RNAs with {args.workers} workers")
    print()
    
    # Process
    rna_args = [(rna_id, seq, output_dir) for rna_id, seq in rnas.items()]
    
    stats = {"msa": 0, "structure": 0}
    start_time = time.time()
    processed = 0
    
    with mp.Pool(processes=args.workers) as pool:
        for result in pool.imap_unordered(process_rna, rna_args):
            processed += 1
            
            if result["msa_success"]:
                stats["msa"] += 1
            if result["structure_success"]:
                stats["structure"] += 1
                completed.add(result["rna_id"])
            
            # Save progress every 10
            if processed % 10 == 0:
                with open(progress_file, 'w') as f:
                    json.dump({
                        "completed": list(completed),
                        "last_updated": datetime.now().isoformat()
                    }, f)
            
            # Print progress
            elapsed = time.time() - start_time
            rate = processed / elapsed if elapsed > 0 else 0
            eta_hours = (len(rnas) - processed) / rate / 3600 if rate > 0 else 0
            
            status = "✓" if result["structure_success"] else "✗"
            method = f"({result['method']})" if result.get("method") else ""
            print(f"[{processed}/{len(rnas)}] {status} {result['rna_id']} {method} | "
                  f"Rate: {rate:.2f}/s | ETA: {eta_hours:.1f}h")
    
    # Final save
    with open(progress_file, 'w') as f:
        json.dump({
            "completed": list(completed),
            "last_updated": datetime.now().isoformat()
        }, f)
    
    # Summary
    total_time = time.time() - start_time
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"RNAs processed: {len(rnas)}")
    print(f"MSAs generated: {stats['msa']}")
    print(f"Structures predicted: {stats['structure']}")
    print(f"Success rate: {stats['structure']/len(rnas)*100:.1f}%")
    print(f"Total time: {total_time/3600:.2f} hours")
    print()


if __name__ == "__main__":
    main()
