import subprocess
from pathlib import Path
import time
import re
import csv

rfam_cm = Path('RNA-Protein/data/rfam/Rfam.cm')
rfam_seed = Path('RNA-Protein/data/rfam/Rfam.seed')


rna_types = {
    'rRNA': [],
    'tRNA': [],
    'miRNA': [],
    'snoRNA': [],
    'snRNA': [],
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
                    rna_type = None
                    
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
                    elif 'lncrna' not in header_lower and 'long non-coding' not in header_lower:
                        # Only use as 'other' if it's NOT a lncRNA
                        rna_type = 'other'
                    
                    if rna_type and len(rna_types[rna_type]) < 5:
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
output_dir = Path('RNA-Protein/data/test_infernal_fair')
output_dir.mkdir(exist_ok=True, parents=True)

# Test each RNA
results = []
success_count = 0

import sys
sys.path.insert(0, 'RNA-Protein/src')

for i, rna in enumerate(test_rnas, 1):
    print(f"[{i}/{len(test_rnas)}] {rna['id']} ({rna['length']} nt, {rna['type']})")
    print(f"  Running cmscan...", end='', flush=True)
    
    # Write query file
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    # Run Infernal cmscan
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
    
    try:
        start = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        elapsed = time.time() - start
        
        # Parse results
        families_found = []
        with open(infernal_out, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 16:
                    family_name = parts[2]
                    family_acc = parts[1]
                    evalue = float(parts[15])
                    families_found.append((family_name, family_acc, evalue))
        
        if families_found:
            best_family = families_found[0][0]
            best_acc = families_found[0][1]
            
            print(f" ✓ Found {best_family} ({best_acc})")
            print(f"    Extracting seed alignment...", end='', flush=True)
            
            msa_generated = False
            a3m_file = None
            n_seqs = 0
            
            try:
                # Extract Stockholm alignment from Rfam.seed
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
                                if re.match(r'^[A-Z0-9]', line):
                                    n_seqs += 1
                                if line.startswith('//'):
                                    break
                
                print(f" {n_seqs} seqs", end='', flush=True)
                
                if n_seqs > 0:
                    print(f", converting to A3M...", end='', flush=True)
                    from rna_msa.alignment_formats import stockholm_to_a3m, validate_a3m
                    
                    a3m_file = output_dir / f"{rna['id']}_{best_acc}.a3m"
                    stockholm_to_a3m(sto_file, a3m_file, query_id=None)
                    
                    stats = validate_a3m(a3m_file)
                    if stats['valid']:
                        msa_generated = True
                        success_count += 1
                        print(f" A3M generated ({stats['n_sequences']} seqs)")
                    else:
                        print(f" A3M invalid")
                else:
                    print(f" Empty seed")
                    
            except Exception as e:
                print(f" Failed: {str(e)[:30]}")
            
            print(f"    Total time: {elapsed:.1f}s")
            
            results.append({
                'rna_id': rna['id'],
                'rna_type': rna['type'],
                'length': rna['length'],
                'success': msa_generated,
                'family': best_family,
                'family_acc': best_acc,
                'n_families': len(families_found),
                'msa_depth': n_seqs,
                'a3m_file': str(a3m_file) if a3m_file else None,
                'time': elapsed
            })
        else:
            print(f" No Rfam family found")
            print(f"    Time: {elapsed:.1f}s")
            
            results.append({
                'rna_id': rna['id'],
                'rna_type': rna['type'],
                'length': rna['length'],
                'success': False,
                'family': None,
                'family_acc': None,
                'n_families': 0,
                'msa_depth': 0,
                'a3m_file': None,
                'time': elapsed
            })
            
    except subprocess.TimeoutExpired:
        print(f" Timeout")
        results.append({
            'rna_id': rna['id'],
            'rna_type': rna['type'],
            'length': rna['length'],
            'success': False,
            'family': None,
            'family_acc': None,
            'n_families': 0,
            'msa_depth': 0,
            'a3m_file': None,
            'time': 120
        })
    except Exception as e:
        print(f" Error: {str(e)[:50]}")
        results.append({
            'rna_id': rna['id'],
            'rna_type': rna['type'],
            'length': rna['length'],
            'success': False,
            'family': None,
            'family_acc': None,
            'n_families': 0,
            'msa_depth': 0,
            'a3m_file': None,
            'time': 0
        })
    
    print()

# Save results
print("=" * 80)
print("RESULTS SUMMARY")
print("=" * 80)
print()

print(f"Total RNAs tested: {len(test_rnas)}")
print(f"Successful MSAs: {success_count}")
print(f"Success rate: {success_count}/{len(test_rnas)} = {success_count/len(test_rnas)*100:.1f}%")
print()

print("Success rate by RNA type:")
for rna_type in ['rRNA', 'tRNA', 'miRNA', 'snoRNA', 'snRNA', 'other']:
    type_results = [r for r in results if r['rna_type'] == rna_type]
    type_success = sum(1 for r in type_results if r['success'])
    if type_results:
        print(f"  {rna_type:8s}: {type_success}/{len(type_results)} ({type_success/len(type_results)*100:.1f}%)")

print()

# Save to TSV
tsv_file = output_dir / 'infernal_fair_results.tsv'
with open(tsv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['rna_id', 'rna_type', 'length', 'success', 'family', 'family_acc', 'n_families', 'msa_depth', 'a3m_file', 'time'], delimiter='\t')
    writer.writeheader()
    writer.writerows(results)

print(f"Results saved to: {tsv_file}")
print()

print("COMPARISON:")
print(f"  Original test (with lncRNAs): 62.9% (22/35)")
print(f"  Fair test (without lncRNAs): {success_count/len(test_rnas)*100:.1f}% ({success_count}/{len(test_rnas)})")
print()

if success_count/len(test_rnas) > 0.7:
    print("Infernal+Rfam performs well on RNA types it's designed for!")
else:
    print("Some RNA types still challenging for Infernal+Rfam")
