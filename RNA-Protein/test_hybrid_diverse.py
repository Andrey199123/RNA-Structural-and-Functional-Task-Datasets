import subprocess
from pathlib import Path
import time
import sys
import csv
import re


# Check databases
rfam_cm = Path('RNA-Protein/data/rfam/Rfam.cm')
rfam_seed = Path('RNA-Protein/data/rfam/Rfam.seed')
rnacentral_db = Path('RNA-Protein/data/rnacentral/rnacentral_active')

rna_types = {
    'rRNA': [],
    'tRNA': [],
    'miRNA': [],
    'snoRNA': [],
    'snRNA': [],
    'lncRNA': [],
    'other': []
}

current_header = None
current_seq = []

with open('RNA-Protein/data/rnacentral/rnacentral_active.fasta', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_header and current_seq:
                seq = ''.join(current_seq)
                seq_len = len(seq)
                
                if 100 <= seq_len <= 2000:
                    seq_id = current_header.split()[0][1:]
                    header_lower = current_header.lower()
                    rna_type = 'other'
                    
                    if 'rrna' in header_lower or 'ribosomal' in header_lower:
                        rna_type = 'rRNA'
                    elif 'trna' in header_lower or 'transfer' in header_lower:
                        rna_type = 'tRNA'
                    elif 'mirna' in header_lower or 'microrna' in header_lower:
                        rna_type = 'miRNA'
                    elif 'snorna' in header_lower:
                        rna_type = 'snoRNA'
                    elif 'snrna' in header_lower:
                        rna_type = 'snRNA'
                    elif 'lncrna' in header_lower or 'long non-coding' in header_lower:
                        rna_type = 'lncRNA'
                    
                    if len(rna_types[rna_type]) < 5:
                        rna_types[rna_type].append({
                            'id': seq_id,
                            'seq': seq,
                            'length': seq_len,
                            'header': current_header,
                            'type': rna_type
                        })
            
            current_header = line
            current_seq = []
            
            if all(len(rnas) >= 5 for rnas in rna_types.values()):
                break
        else:
            current_seq.append(line)

test_rnas = []
for rna_type, rnas in rna_types.items():
    test_rnas.extend(rnas)
    print(f"  {rna_type}: {len(rnas)} RNAs")
  
# Create output directory
output_dir = Path('RNA-Protein/data/test_hybrid_diverse')
output_dir.mkdir(exist_ok=True, parents=True)

# Test each RNA with hybrid approach
results = []
infernal_success = 0
blast_success = 0
hybrid_success = 0

for i, rna in enumerate(test_rnas, 1):
    print(f"[{i}/{len(test_rnas)}] {rna['id']} ({rna['length']} nt, {rna['type']})")
    
    start_time = time.time()
    
    # STEP 1: Try Infernal+Rfam
    print(f"  [1] Trying Infernal+Rfam...", end='', flush=True)
    
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    infernal_success_flag = False
    blast_success_flag = False
    final_method = None
    a3m_file = None
    
    try:
        # Run cmscan
        infernal_out = output_dir / f"{rna['id']}_infernal.tblout"
        cmd = [
            '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/cmscan',
            '--tblout', str(infernal_out),
            '--cut_ga',
            '--nohmmonly',
            '--cpu', '4',
            str(rfam_cm),
            str(query_file)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
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
        
        if families_found:
            best_acc = families_found[0]
            print(f" found {best_acc}", end='', flush=True)
            
            # Extract seed alignment
            sys.path.insert(0, 'RNA-Protein/src')
            from rna_msa.alignment_formats import stockholm_to_a3m, validate_a3m
            
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
            
            # Convert to A3M
            a3m_file = output_dir / f"{rna['id']}_{best_acc}.a3m"
            stockholm_to_a3m(sto_file, a3m_file, query_id=None)
            
            stats = validate_a3m(a3m_file)
            if stats['valid']:
                infernal_success_flag = True
                infernal_success += 1
                hybrid_success += 1
                final_method = 'Infernal'
                print(f" Success ({stats['n_sequences']} seqs)")
            else:
                print(f" Invalid A3M", end='')
        else:
            print(f" no family", end='')
    
    except Exception as e:
        print(f" Error", end='')
    
    # STEP 2: If Infernal failed, try BLAST
    if not infernal_success_flag:
        print(f"\n  [2] Trying BLAST+RNAcentral...", end='', flush=True)
        
        try:
            # Run BLAST
            blast_out = output_dir / f"{rna['id']}_blast.txt"
            cmd = [
                '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/blastn',
                '-query', str(query_file),
                '-db', 'RNA-Protein/data/rnacentral/rnacentral_active',
                '-evalue', '0.001',
                '-outfmt', '6 sseqid pident length',
                '-max_target_seqs', '100',
                '-num_threads', '4'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            homolog_ids = set()
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        seq_id = parts[0]
                        pident = float(parts[1])
                        if pident >= 70:
                            homolog_ids.add(seq_id)
            
            n_homologs = len(homolog_ids)
            print(f" {n_homologs} homologs", end='', flush=True)
            
            if n_homologs >= 10:
                # Extract and align
                id_file = output_dir / f"{rna['id']}_ids.txt"
                with open(id_file, 'w') as f:
                    for seq_id in homolog_ids:
                        f.write(f"{seq_id}\n")
                
                homologs_file = output_dir / f"{rna['id']}_homologs.fa"
                cmd = [
                    '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/blastdbcmd',
                    '-db', 'RNA-Protein/data/rnacentral/rnacentral_active',
                    '-entry_batch', str(id_file),
                    '-line_length', '100000',
                    '-out', str(homologs_file)
                ]
                subprocess.run(cmd, capture_output=True, text=True, timeout=60)
                
                aligned_file = output_dir / f"{rna['id']}_aligned.fa"
                cmd = ['/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/mafft', '--auto', '--thread', '4', str(homologs_file)]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                
                with open(aligned_file, 'w') as f:
                    f.write(result.stdout)
                
                from rna_msa.alignment_formats import fasta_to_a3m, validate_a3m
                
                a3m_file = output_dir / f"{rna['id']}_blast.a3m"
                fasta_to_a3m(aligned_file, a3m_file, query_id=rna['id'])
                
                stats = validate_a3m(a3m_file)
                if stats['valid']:
                    blast_success_flag = True
                    blast_success += 1
                    hybrid_success += 1
                    final_method = 'BLAST'
                    print(f" ✓ Success ({stats['n_sequences']} seqs)")
                else:
                    print(f" ✗ Invalid A3M")
            else:
                print(f" ✗ Insufficient")
        
        except Exception as e:
            print(f" ✗ Error")
    
    elapsed = time.time() - start_time
    
    if not infernal_success_flag and not blast_success_flag:
        print(f"  [FINAL] ✗ Both methods failed")
    
    print(f"  Time: {elapsed:.1f}s")
    print()
    
    results.append({
        'rna_id': rna['id'],
        'rna_type': rna['type'],
        'length': rna['length'],
        'infernal_success': infernal_success_flag,
        'blast_success': blast_success_flag,
        'hybrid_success': infernal_success_flag or blast_success_flag,
        'final_method': final_method,
        'a3m_file': str(a3m_file) if a3m_file else None,
        'time': elapsed
    })


print(f"Total RNAs tested: {len(test_rnas)}")
print(f"Infernal-only success: {infernal_success}/{len(test_rnas)} ({infernal_success/len(test_rnas)*100:.1f}%)")
print(f"BLAST fallback used: {blast_success}/{len(test_rnas)} ({blast_success/len(test_rnas)*100:.1f}%)")
print(f"HYBRID total success: {hybrid_success}/{len(test_rnas)} ({hybrid_success/len(test_rnas)*100:.1f}%)")

print("Success rate by RNA type:")
for rna_type in ['rRNA', 'tRNA', 'miRNA', 'snoRNA', 'snRNA', 'lncRNA', 'other']:
    type_results = [r for r in results if r['rna_type'] == rna_type]
    type_success = sum(1 for r in type_results if r['hybrid_success'])
    if type_results:
        print(f"  {rna_type:8s}: {type_success}/{len(type_results)} ({type_success/len(type_results)*100:.1f}%)")

print()

# Save to TSV
tsv_file = output_dir / 'hybrid_diverse_results.tsv'
with open(tsv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['rna_id', 'rna_type', 'length', 'infernal_success', 'blast_success', 'hybrid_success', 'final_method', 'a3m_file', 'time'], delimiter='\t')
    writer.writeheader()
    writer.writerows(results)

print(f"Results saved to: {tsv_file}")
print("COMPARISON:")
print(f"  Infernal-only: 62.9% (from previous test)")
print(f"  BLAST-only: 74.3% (from previous test)")
print(f"  HYBRID: {hybrid_success/len(test_rnas)*100:.1f}% (this test)")


