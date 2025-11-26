#!/usr/bin/env python3
"""
Test Infernal + Rfam on the same 30 diverse RNAs from benchmark
Compare results with existing BLAST benchmark
"""

import subprocess
from pathlib import Path
import time

print("=" * 80)
print("TESTING INFERNAL + RFAM")
print("On 30 Diverse RNAs (same as BLAST benchmark)")
print("=" * 80)
print()

# Use the same 30 RNAs from benchmark
benchmark_dir = Path('data/msa_benchmark')
query_files = sorted(benchmark_dir.glob('*_query.fa'))[:30]

print(f"Found {len(query_files)} test RNAs")
print()

# Check Rfam database
rfam_cm = Path('data/rfam/Rfam.cm')
rfam_seed = Path('data/rfam/Rfam.seed')

if not rfam_cm.exists():
    print(f"✗ Rfam.cm not found at {rfam_cm}")
    exit(1)

if not rfam_seed.exists():
    print(f"✗ Rfam.seed not found at {rfam_seed}")
    exit(1)

print("✓ Rfam database found")
print()

# Results storage
results = []

print("Testing each RNA with Infernal + Rfam...")
print()

for i, query_file in enumerate(query_files, 1):
    rna_id = query_file.stem.replace('_query', '')
    
    # Read sequence to get length
    with open(query_file, 'r') as f:
        lines = f.readlines()
        seq = ''.join(line.strip() for line in lines[1:])
        seq_len = len(seq)
    
    print(f"[{i}/30] {rna_id} ({seq_len} nt)")
    print(f"  Running cmscan...", end='', flush=True)
    
    # Run Infernal cmscan
    output_dir = Path('data/test_infernal_rfam')
    output_dir.mkdir(exist_ok=True)
    
    infernal_out = output_dir / f'{rna_id}_infernal.tblout'
    
    cmd = [
        'cmscan',
        '--tblout', str(infernal_out),
        '--cut_ga',  # Use gathering threshold
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
                    family_name = parts[2]  # Rfam ID
                    family_acc = parts[1]   # Rfam accession
                    evalue = float(parts[15])
                    families_found.append((family_name, family_acc, evalue))
        
        if families_found:
            best_family = families_found[0][0]
            best_acc = families_found[0][1]
            
            print(f" ✓ Found {len(families_found)} family match(es)")
            print(f"    Best: {best_family} ({best_acc})")
            
            # Get seed alignment size for best family
            print(f"    Counting sequences in seed alignment...", end='', flush=True)
            try:
                # Count sequence lines using Python (handles encoding issues)
                import re
                n_seqs = 0
                in_family = False
                with open(rfam_seed, 'r', errors='ignore') as f:
                    for line in f:
                        if line.startswith(f'#=GF AC   {best_acc}'):
                            in_family = True
                        elif line.startswith('//') and in_family:
                            break
                        elif in_family and re.match(r'^[A-Z0-9]', line):
                            n_seqs += 1
                
                if n_seqs > 0:
                    print(f" {n_seqs} sequences")
                else:
                    print(f" 0 sequences (empty seed)")
            except Exception as e:
                print(f" Error: {str(e)[:30]}")
                n_seqs = 0
            
            results.append({
                'rna_id': rna_id,
                'length': seq_len,
                'success': True,
                'family': best_family,
                'family_acc': best_acc,
                'n_families': len(families_found),
                'msa_depth': n_seqs,
                'time': elapsed
            })
            
            print(f"    Time: {elapsed:.1f}s")
        else:
            print(f" ✗ No Rfam family found")
            print(f"    Time: {elapsed:.1f}s")
            
            results.append({
                'rna_id': rna_id,
                'length': seq_len,
                'success': False,
                'family': None,
                'family_acc': None,
                'n_families': 0,
                'msa_depth': 0,
                'time': elapsed
            })
            
    except subprocess.TimeoutExpired:
        print(f" ✗ Timeout (>120s)")
        results.append({
            'rna_id': rna_id,
            'length': seq_len,
            'success': False,
            'family': None,
            'family_acc': None,
            'n_families': 0,
            'msa_depth': 0,
            'time': 120
        })
    except Exception as e:
        print(f" ✗ Error: {str(e)[:50]}")
        results.append({
            'rna_id': rna_id,
            'length': seq_len,
            'success': False,
            'family': None,
            'family_acc': None,
            'n_families': 0,
            'msa_depth': 0,
            'time': 0
        })
    
    print()

# ========================================
# Summary Statistics
# ========================================
print("=" * 80)
print("INFERNAL + RFAM RESULTS")
print("=" * 80)
print()

total = len(results)
successful = [r for r in results if r['success']]
n_success = len(successful)
success_rate = n_success / total * 100

print(f"Success Rate: {n_success}/{total} ({success_rate:.1f}%)")
print(f"  (Success = found Rfam family with seed alignment)")
print()

if n_success > 0:
    # MSA depth statistics
    depths = [r['msa_depth'] for r in successful if r['msa_depth'] > 0]
    if depths:
        print(f"MSA Depth (seed alignment size):")
        print(f"  Min:     {min(depths)} sequences")
        print(f"  Max:     {max(depths)} sequences")
        print(f"  Average: {sum(depths)/len(depths):.1f} sequences")
        print(f"  Median:  {sorted(depths)[len(depths)//2]} sequences")
        print()
    
    # Time statistics
    times = [r['time'] for r in successful]
    print(f"Search Time:")
    print(f"  Min:     {min(times):.1f}s")
    print(f"  Max:     {max(times):.1f}s")
    print(f"  Average: {sum(times)/len(times):.1f}s")
    print()
    
    # Families found
    families = {}
    for r in successful:
        fam = r['family']
        if fam not in families:
            families[fam] = []
        families[fam].append(r['rna_id'])
    
    print(f"Rfam Families Found: {len(families)}")
    for fam, rnas in sorted(families.items(), key=lambda x: -len(x[1])):
        print(f"  {fam}: {len(rnas)} RNAs")
        if len(rnas) <= 3:
            print(f"    ({', '.join(rnas)})")
    print()

# Failed RNAs
failed = [r for r in results if not r['success']]
if failed:
    print(f"Failed RNAs ({len(failed)}):")
    for r in failed:
        print(f"  {r['rna_id']} ({r['length']} nt) - No Rfam family found")
    print()

# ========================================
# Compare with BLAST benchmark
# ========================================
print("=" * 80)
print("COMPARISON WITH BLAST BENCHMARK")
print("=" * 80)
print()

# Load BLAST benchmark results
blast_file = Path('data/msa_benchmark/benchmark_results.tsv')
if blast_file.exists():
    # Parse TSV manually
    blast_data = {}
    with open(blast_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                rna_id = parts[0]
                homologs = int(parts[3])
                success = parts[7] == 'True'
                blast_data[rna_id] = {'homologs': homologs, 'success': success}
    
    # Match by RNA_ID
    blast_success = 0
    blast_total = 0
    blast_depths = []
    
    for r in results:
        if r['rna_id'] in blast_data:
            blast_total += 1
            if blast_data[r['rna_id']]['success']:
                blast_success += 1
                blast_depths.append(blast_data[r['rna_id']]['homologs'])
    
    blast_rate = blast_success / blast_total * 100 if blast_total > 0 else 0
    
    print("Success Rate:")
    print(f"  Infernal+Rfam:     {success_rate:5.1f}% ({n_success}/{total})")
    print(f"  BLAST+RNAcentral:  {blast_rate:5.1f}% ({blast_success}/{blast_total})")
    print()
    
    if depths and blast_depths:
        print("Average MSA Depth:")
        print(f"  Infernal+Rfam:     {sum(depths)/len(depths):5.1f} sequences")
        print(f"  BLAST+RNAcentral:  {sum(blast_depths)/len(blast_depths):5.1f} sequences")
        print()
else:
    print("⚠ BLAST benchmark results not found")
    print()

# ========================================
# Save results
# ========================================
results_file = Path('data/test_infernal_rfam/infernal_results.tsv')
with open(results_file, 'w') as f:
    f.write("RNA_ID\tLength\tSuccess\tFamily\tFamily_Acc\tN_Families\tMSA_Depth\tTime\n")
    for r in results:
        f.write(f"{r['rna_id']}\t{r['length']}\t{r['success']}\t{r.get('family', 'N/A')}\t{r.get('family_acc', 'N/A')}\t{r['n_families']}\t{r['msa_depth']}\t{r['time']:.1f}\n")

print(f"Detailed results saved to: {results_file}")
print()

# ========================================
# Conclusion
# ========================================
print("=" * 80)
print("CONCLUSION")
print("=" * 80)
print()

if success_rate >= 90:
    print("✓ EXCELLENT: Infernal+Rfam works very well for these RNAs")
    print(f"  {success_rate:.0f}% of RNAs found in Rfam families")
elif success_rate >= 70:
    print("✓ GOOD: Infernal+Rfam works for most RNAs")
    print(f"  {success_rate:.0f}% of RNAs found in Rfam families")
elif success_rate >= 50:
    print("⚠ FAIR: Infernal+Rfam works for some RNAs")
    print(f"  {success_rate:.0f}% of RNAs found in Rfam families")
else:
    print("✗ LIMITED: Infernal+Rfam works for few RNAs")
    print(f"  Only {success_rate:.0f}% of RNAs found in Rfam families")

print()
print("Key findings:")
print(f"  • {success_rate:.0f}% of test RNAs are in Rfam families")
if depths:
    print(f"  • Average seed alignment: {sum(depths)/len(depths):.0f} sequences")
if times:
    print(f"  • Average search time: {sum(times)/len(times):.0f}s")
print()
