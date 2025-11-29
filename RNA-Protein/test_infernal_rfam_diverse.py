#!/usr/bin/env python3
"""
Tests Infernal + Rfam on diverse RNA types from RNAcentral.
"""

import subprocess
from pathlib import Path
import time
import re

print("TESTING INFERNAL + RFAM ON DIVERSE RNA TYPES")
print()

rfam_cm = Path('data/rfam/Rfam.cm')
rfam_seed = Path('data/rfam/Rfam.seed')

if not rfam_cm.exists() or not rfam_seed.exists():
    print("Rfam database not found")
    exit(1)

print("Rfam database found")
print()

print("Extracting diverse RNA sample from RNAcentral...")
print("Looking for different RNA types: rRNA, tRNA, miRNA, snoRNA, lncRNA, etc.")
print()

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

with open('data/rnacentral/rnacentral_active.fasta', 'r') as f:
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
                            'header': current_header
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

print(f"\nTotal test RNAs: {len(test_rnas)}")
print()

# Create output directory
output_dir = Path('data/test_infernal_diverse')
output_dir.mkdir(exist_ok=True)

# Test each RNA
results = []

print("Testing each RNA with Infernal + Rfam...")
print()

for i, rna in enumerate(test_rnas, 1):
    print(f"[{i}/{len(test_rnas)}] {rna['id']} ({rna['length']} nt)")
    
    # Determine RNA type from header
    header_lower = rna['header'].lower()
    if 'rrna' in header_lower:
        rna_type = 'rRNA'
    elif 'trna' in header_lower:
        rna_type = 'tRNA'
    elif 'mirna' in header_lower:
        rna_type = 'miRNA'
    elif 'snorna' in header_lower:
        rna_type = 'snoRNA'
    elif 'snrna' in header_lower:
        rna_type = 'snRNA'
    elif 'lncrna' in header_lower:
        rna_type = 'lncRNA'
    else:
        rna_type = 'other'
    
    print(f"  Type: {rna_type}")
    print(f"  Running cmscan...", end='', flush=True)
    
    # Write query file
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    # Run Infernal cmscan
    infernal_out = output_dir / f"{rna['id']}_infernal.tblout"
    
    cmd = [
        'cmscan',
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
            
            print(f"Found {len(families_found)} family match(es)")
            print(f"    Best: {best_family} ({best_acc})")
            
            # Extract seed alignment and convert to A3M
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
                    # Convert to A3M (use first sequence as query since our query isn't in the seed)
                    print(f", converting to A3M...", end='', flush=True)
                    from src.rna_msa.alignment_formats import stockholm_to_a3m, validate_a3m
                    
                    a3m_file = output_dir / f"{rna['id']}_{best_acc}.a3m"
                    # Use None for query_id to use first sequence in alignment
                    stockholm_to_a3m(sto_file, a3m_file, query_id=None)
                    
                    # Validate
                    stats = validate_a3m(a3m_file)
                    if stats['valid']:
                        msa_generated = True
                        print(f"A3M generated ({stats['n_sequences']} seqs)")
                    else:
                        print(f"A3M invalid")
                else:
                    print(f"Empty seed")
                    
            except Exception as e:
                print(f"Failed: {str(e)[:30]}")
            
            print(f"    Total time: {elapsed:.1f}s")
            
            results.append({
                'rna_id': rna['id'],
                'rna_type': rna_type,
                'length': rna['length'],
                'success': msa_generated,
                'msa_generated': msa_generated,
                'family': best_family,
                'family_acc': best_acc,
                'n_families': len(families_found),
                'msa_depth': n_seqs,
                'a3m_file': str(a3m_file) if a3m_file else None,
                'time': elapsed
            })
        else:
            print(f"No Rfam family found")
            print(f"    Time: {elapsed:.1f}s")
            
            results.append({
                'rna_id': rna['id'],
                'rna_type': rna_type,
                'length': rna['length'],
                'success': False,
                'msa_generated': False,
                'family': None,
                'family_acc': None,
                'n_families': 0,
                'msa_depth': 0,
                'a3m_file': None,
                'time': elapsed
            })
            
    except subprocess.TimeoutExpired:
        print(f"Timeout")
        results.append({
            'rna_id': rna['id'],
            'rna_type': rna_type,
            'length': rna['length'],
            'success': False,
            'msa_generated': False,
            'family': None,
            'family_acc': None,
            'n_families': 0,
            'msa_depth': 0,
            'a3m_file': None,
            'time': 120
        })
    except Exception as e:
        print(f"Error: {str(e)[:50]}")
        results.append({
            'rna_id': rna['id'],
            'rna_type': rna_type,
            'length': rna['length'],
            'success': False,
            'msa_generated': False,
            'family': None,
            'family_acc': None,
            'n_families': 0,
            'msa_depth': 0,
            'a3m_file': None,
            'time': 0
        })
    
    print()

print("RESULTS BY RNA TYPE")
print()

# Group by RNA type
by_type = {}
for r in results:
    rna_type = r['rna_type']
    if rna_type not in by_type:
        by_type[rna_type] = []
    by_type[rna_type].append(r)

for rna_type in sorted(by_type.keys()):
    type_results = by_type[rna_type]
    total = len(type_results)
    successful = [r for r in type_results if r['success']]
    n_success = len(successful)
    success_rate = n_success / total * 100 if total > 0 else 0
    
    print(f"{rna_type}:")
    print(f"  Tested: {total}")
    print(f"  Success: {n_success}/{total} ({success_rate:.1f}%)")
    
    if successful:
        depths = [r['msa_depth'] for r in successful if r['msa_depth'] > 0]
        if depths:
            print(f"  Avg MSA depth: {sum(depths)/len(depths):.1f} sequences")
        
        families = set(r['family'] for r in successful)
        print(f"  Families found: {', '.join(sorted(families))}")
    
    print()

print("OVERALL SUMMARY")
print()

total = len(results)
msa_generated = sum(1 for r in results if r.get('msa_generated', False))
msa_rate = msa_generated / total * 100

print(f"Total RNAs tested: {total}")
print(f"MSAs generated: {msa_generated}/{total} ({msa_rate:.1f}%)")
print(f"  (Success = A3M file generated)")
print()

if successful:
    depths = [r['msa_depth'] for r in successful if r['msa_depth'] > 0]
    if depths:
        print(f"MSA Depth:")
        print(f"  Min: {min(depths)} sequences")
        print(f"  Max: {max(depths)} sequences")
        print(f"  Avg: {sum(depths)/len(depths):.1f} sequences")
        print()
    
    times = [r['time'] for r in successful]
    print(f"Search Time:")
    print(f"  Avg: {sum(times)/len(times):.1f}s")
    print()

# Save results
results_file = output_dir / 'diverse_test_results.tsv'
with open(results_file, 'w') as f:
    f.write("RNA_ID\tRNA_Type\tLength\tMSA_Generated\tFamily\tFamily_Acc\tN_Families\tMSA_Depth\tA3M_File\tTime\n")
    for r in results:
        f.write(f"{r['rna_id']}\t{r['rna_type']}\t{r['length']}\t{r.get('msa_generated', False)}\t{r.get('family', 'N/A')}\t{r.get('family_acc', 'N/A')}\t{r['n_families']}\t{r['msa_depth']}\t{r.get('a3m_file', 'N/A')}\t{r['time']:.1f}\n")

print(f"Detailed results saved to: {results_file}")
