#!/usr/bin/env python3
"""Process 230 Rfam sequences with 12 threads"""

import subprocess
from pathlib import Path
import time
import sys
import csv

sys.path.insert(0, str(Path(__file__).parent / "RNA-Protein/src"))
from rna_msa.alignment_formats import fasta_to_a3m, stockholm_to_a3m, validate_a3m


def parse_fasta(fasta_file):
    rnas = []
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id and current_seq:
                    seq = ''.join(current_seq)
                    rnas.append({'id': current_id, 'seq': seq, 'length': len(seq)})
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id and current_seq:
            seq = ''.join(current_seq)
            rnas.append({'id': current_id, 'seq': seq, 'length': len(seq)})
    
    return rnas


def run_infernal_rfam(rna, output_dir, rfam_cm, rfam_seed, cmscan_path, threads=12):
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    try:
        infernal_out = output_dir / f"{rna['id']}_infernal.tblout"
        cmd = [cmscan_path, '--tblout', str(infernal_out), '--cut_ga', '--nohmmonly',
               '--cpu', str(threads), str(rfam_cm), str(query_file)]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        families_found = []
        with open(infernal_out, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 16:
                        families_found.append(parts[1])
        
        if not families_found:
            return False, None, 0, "No Rfam family found"
        
        best_acc = families_found[0]
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
        
        a3m_file = output_dir / f"{rna['id']}_{best_acc}.a3m"
        stockholm_to_a3m(sto_file, a3m_file, query_id=rna['id'], query_sequence=rna['seq'])
        
        stats = validate_a3m(a3m_file)
        if stats['valid']:
            return True, a3m_file, stats['n_sequences'], None
        else:
            return False, None, 0, f"Invalid A3M: {stats.get('error', 'Unknown')}"
    
    except Exception as e:
        return False, None, 0, str(e)


def run_blast_mafft(rna, output_dir, rnacentral_db, blastn_path, blastdbcmd_path, 
                    mafft_path, threads=12):
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    try:
        blast_out = output_dir / f"{rna['id']}_blast.txt"
        cmd = [blastn_path, '-query', str(query_file), '-db', str(rnacentral_db),
               '-evalue', '0.001', '-outfmt', '6 sseqid pident length',
               '-max_target_seqs', '100', '-num_threads', str(threads)]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        homolog_ids = set()
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 3 and float(parts[1]) >= 70.0:
                    homolog_ids.add(parts[0])
        
        if len(homolog_ids) < 10:
            return False, None, 0, f"Insufficient homologs ({len(homolog_ids)} < 10)"
        
        id_file = output_dir / f"{rna['id']}_ids.txt"
        with open(id_file, 'w') as f:
            for seq_id in homolog_ids:
                f.write(f"{seq_id}\n")
        
        homologs_file = output_dir / f"{rna['id']}_homologs.fa"
        cmd = [blastdbcmd_path, '-db', str(rnacentral_db), '-entry_batch', str(id_file),
               '-line_length', '100000', '-out', str(homologs_file)]
        subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        with open(homologs_file, 'r') as f:
            homolog_content = f.read()
        
        with open(homologs_file, 'w') as f:
            f.write(f">{rna['id']}\n{rna['seq']}\n")
            f.write(homolog_content)
        
        aligned_file = output_dir / f"{rna['id']}_aligned.fa"
        cmd = [mafft_path, '--auto', '--thread', str(threads), str(homologs_file)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        
        with open(aligned_file, 'w') as f:
            f.write(result.stdout)
        
        a3m_file = output_dir / f"{rna['id']}_blast.a3m"
        fasta_to_a3m(aligned_file, a3m_file, query_id=rna['id'])
        
        stats = validate_a3m(a3m_file)
        if stats['valid']:
            return True, a3m_file, stats['n_sequences'], None
        else:
            return False, None, 0, f"Invalid A3M: {stats.get('error', 'Unknown')}"
    
    except Exception as e:
        return False, None, 0, str(e)


# Main processing
input_fasta = Path('sent/Rfam_sequences.txt')
output_dir = Path('RNA-Protein/data/rfam_230_output')
output_dir.mkdir(exist_ok=True, parents=True)

rfam_cm = Path('RNA-Protein/data/rfam/Rfam.cm')
rfam_seed = Path('RNA-Protein/data/rfam/Rfam.seed')
rnacentral_db = Path('RNA-Protein/data/rnacentral/rnacentral_active')

cmscan_path = '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/cmscan'
blastn_path = '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/blastn'
blastdbcmd_path = '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/blastdbcmd'
mafft_path = '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/mafft'

print("Loading RNA sequences...")
rnas = parse_fasta(input_fasta)
print(f"Loaded {len(rnas)} RNA sequences\n")

results = []
infernal_success = 0
blast_success = 0
total_success = 0

start_time_total = time.time()

for i, rna in enumerate(rnas, 1):
    print(f"[{i}/{len(rnas)}] {rna['id']} ({rna['length']} nt)")
    
    start_time = time.time()
    
    print(f"  [1] Trying Infernal+Rfam...", end='', flush=True)
    success, a3m_file, n_seqs, error = run_infernal_rfam(
        rna, output_dir, rfam_cm, rfam_seed, cmscan_path, threads=12
    )
    
    if success:
        infernal_success += 1
        total_success += 1
        final_method = 'Infernal'
        print(f" ✓ Success ({n_seqs} seqs)")
    else:
        print(f" ✗ {error}")
        
        print(f"  [2] Trying BLAST+MAFFT...", end='', flush=True)
        success, a3m_file, n_seqs, error = run_blast_mafft(
            rna, output_dir, rnacentral_db, blastn_path, 
            blastdbcmd_path, mafft_path, threads=12
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
    
    print(f"  Time: {elapsed:.1f}s\n")
    
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

elapsed_total = time.time() - start_time_total

print(f"\n{'='*70}")
print(f"FINAL SUMMARY")
print(f"{'='*70}")
print(f"Total RNAs processed: {len(rnas)}")
print(f"Infernal success: {infernal_success}/{len(rnas)} ({infernal_success/len(rnas)*100:.1f}%)")
print(f"BLAST success: {blast_success}/{len(rnas)} ({blast_success/len(rnas)*100:.1f}%)")
print(f"HYBRID total success: {total_success}/{len(rnas)} ({total_success/len(rnas)*100:.1f}%)")
print(f"Failed: {len(rnas) - total_success}/{len(rnas)} ({(len(rnas) - total_success)/len(rnas)*100:.1f}%)")
print(f"Total time: {elapsed_total/60:.1f} minutes")
print(f"Average time per RNA: {elapsed_total/len(rnas):.1f} seconds")

tsv_file = output_dir / 'rfam_230_results.tsv'
with open(tsv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=[
        'rna_id', 'length', 'infernal_success', 'blast_success', 
        'hybrid_success', 'final_method', 'a3m_file', 'n_sequences', 'time'
    ], delimiter='\t')
    writer.writeheader()
    writer.writerows(results)

print(f"\nResults saved to: {tsv_file}")
print(f"A3M files saved to: {output_dir}")
