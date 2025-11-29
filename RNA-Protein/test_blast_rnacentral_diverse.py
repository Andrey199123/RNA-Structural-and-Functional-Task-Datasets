#!/usr/bin/env python3
"""
Tests BLAST + RNAcentral on the same diverse RNAs as the Infernal test.
"""

import subprocess
from pathlib import Path
import time

print("=" * 80)
print("TESTING BLAST + RNACENTRAL ON DIVERSE RNA TYPES")
print("Using the SAME RNAs as Infernal test for fair comparison")
print("=" * 80)
print()

infernal_dir = Path('data/test_infernal_diverse')
if not infernal_dir.exists():
    print("✗ Infernal test results not found")
    print("  Run test_infernal_rfam_diverse.py first")
    exit(1)

# Read the query files from Infernal test
query_files = sorted(infernal_dir.glob('*_query.fa'))

if not query_files:
    print("✗ No query files found")
    exit(1)

print(f"Found {len(query_files)} test RNAs from Infernal test")
print()

infernal_results = {}
with open(infernal_dir / 'diverse_test_results.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            rna_id = parts[0]
            rna_type = parts[1]
            length = int(parts[2])
            infernal_results[rna_id] = {
                'type': rna_type,
                'length': length
            }

print("RNA types in test:")
from collections import Counter
type_counts = Counter(r['type'] for r in infernal_results.values())
for rna_type, count in sorted(type_counts.items()):
    print(f"  {rna_type}: {count}")
print()

output_dir = Path('data/test_blast_diverse')
output_dir.mkdir(exist_ok=True)

results = []

print("Testing each RNA with BLAST + RNAcentral...")
print()

for i, query_file in enumerate(query_files, 1):
    rna_id = query_file.stem.replace('_query', '')
    
    if rna_id not in infernal_results:
        continue
    
    rna_info = infernal_results[rna_id]
    rna_type = rna_info['type']
    length = rna_info['length']
    
    print(f"[{i}/{len(query_files)}] {rna_id} ({length} nt)")
    print(f"  Type: {rna_type}")
    print(f"  Running BLAST...", end='', flush=True)
    
    cmd = [
        'blastn',
        '-query', str(query_file),
        '-db', 'data/rnacentral/rnacentral_active',
        '-evalue', '0.001',
        '-outfmt', '6 sseqid pident length',
        '-max_target_seqs', '100',
        '-num_threads', '4'
    ]
    
    try:
        start = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60)
        elapsed = time.time() - start
        
        hits = []
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 3:
                    seq_id = parts[0]
                    identity = float(parts[1])
                    align_len = int(parts[2])
                    
                    # Quality filter: ≥70% identity, ≥50% of query length
                    if identity >= 70 and align_len >= length * 0.5:
                        hits.append((seq_id, identity, align_len))
        
        n_homologs = len(hits)
        avg_identity = sum(h[1] for h in hits) / len(hits) if hits else 0
        
        # Generate actual MSA if we have enough homologs
        msa_generated = False
        a3m_file = None
        
        if n_homologs >= 5:
            print(f" ✓ {n_homologs} homologs (avg {avg_identity:.1f}% identity)")
            print(f"    Generating MSA...", end='', flush=True)
            
            try:
                # Extract homolog sequences
                id_file = output_dir / f"{rna_id}_ids.txt"
                with open(id_file, 'w') as f:
                    for hit in hits[:50]:  # Limit to 50 for speed
                        f.write(f"{hit[0]}\n")
                
                homologs_file = output_dir / f"{rna_id}_homologs.fa"
                cmd = [
                    'blastdbcmd',
                    '-db', 'data/rnacentral/rnacentral_active',
                    '-entry_batch', str(id_file),
                    '-out', str(homologs_file)
                ]
                subprocess.run(cmd, capture_output=True, check=True, timeout=30)
                
                # Align with MAFFT
                aligned_file = output_dir / f"{rna_id}_aligned.fa"
                cmd = ['mafft', '--auto', '--quiet', '--thread', '4', str(homologs_file)]
                result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=120)
                with open(aligned_file, 'w') as f:
                    f.write(result.stdout)
                
                # Convert to A3M
                from src.rna_msa.alignment_formats import fasta_to_a3m, validate_a3m
                a3m_file = output_dir / f"{rna_id}_aligned.a3m"
                fasta_to_a3m(aligned_file, a3m_file, query_id=rna_id)
                
                # Validate
                stats = validate_a3m(a3m_file)
                if stats['valid']:
                    msa_generated = True
                    print(f" ✓ A3M generated ({stats['n_sequences']} seqs)")
                else:
                    print(f" ✗ A3M invalid")
                    
            except Exception as e:
                print(f" ✗ MSA failed: {str(e)[:30]}")
            
            print(f"    Total time: {elapsed:.1f}s")
            success = msa_generated
        elif n_homologs > 0:
            print(f" ⚠ Only {n_homologs} homologs (need ≥5, {elapsed:.1f}s)")
            success = False
        else:
            print(f" ✗ No homologs found ({elapsed:.1f}s)")
            success = False
        
        results.append({
            'rna_id': rna_id,
            'rna_type': rna_type,
            'length': length,
            'success': success,
            'msa_generated': msa_generated,
            'n_homologs': n_homologs,
            'avg_identity': avg_identity,
            'a3m_file': str(a3m_file) if a3m_file else None,
            'time': elapsed
        })
        
    except subprocess.TimeoutExpired:
        print(f" ✗ Timeout")
        results.append({
            'rna_id': rna_id,
            'rna_type': rna_type,
            'length': length,
            'success': False,
            'n_homologs': 0,
            'avg_identity': 0,
            'time': 60
        })
    except Exception as e:
        print(f" ✗ Error: {str(e)[:50]}")
        results.append({
            'rna_id': rna_id,
            'rna_type': rna_type,
            'length': length,
            'success': False,
            'n_homologs': 0,
            'avg_identity': 0,
            'time': 0
        })
    
    print()

# Summary by RNA type
print("=" * 80)
print("RESULTS BY RNA TYPE")
print("=" * 80)
print()

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
        avg_depth = sum(r['n_homologs'] for r in successful) / len(successful)
        avg_identity = sum(r['avg_identity'] for r in successful) / len(successful)
        print(f"  Avg homologs: {avg_depth:.1f}")
        print(f"  Avg identity: {avg_identity:.1f}%")
    
    print()

# Overall summary
print("=" * 80)
print("OVERALL SUMMARY")
print("=" * 80)
print()

total = len(results)
successful = [r for r in results if r['success']]
n_success = len(successful)
success_rate = n_success / total * 100

msa_generated = sum(1 for r in results if r.get('msa_generated', False))

print(f"Total RNAs tested: {total}")
print(f"MSAs generated: {msa_generated}/{total} ({msa_generated/total*100:.1f}%)")
print(f"  (Success = A3M file generated)")
print()

if successful:
    homolog_counts = [r['n_homologs'] for r in successful]
    identities = [r['avg_identity'] for r in successful]
    times = [r['time'] for r in successful]
    
    print(f"Homolog Count:")
    print(f"  Min: {min(homolog_counts)}")
    print(f"  Max: {max(homolog_counts)}")
    print(f"  Avg: {sum(homolog_counts)/len(homolog_counts):.1f}")
    print()
    
    print(f"Sequence Identity:")
    print(f"  Min: {min(identities):.1f}%")
    print(f"  Max: {max(identities):.1f}%")
    print(f"  Avg: {sum(identities)/len(identities):.1f}%")
    print()
    
    print(f"Search Time:")
    print(f"  Avg: {sum(times)/len(times):.1f}s")
    print()

# Save results
results_file = output_dir / 'diverse_test_results.tsv'
with open(results_file, 'w') as f:
    f.write("RNA_ID\tRNA_Type\tLength\tMSA_Generated\tN_Homologs\tAvg_Identity\tA3M_File\tTime\n")
    for r in results:
        f.write(f"{r['rna_id']}\t{r['rna_type']}\t{r['length']}\t{r.get('msa_generated', False)}\t{r['n_homologs']}\t{r['avg_identity']:.1f}\t{r.get('a3m_file', 'N/A')}\t{r['time']:.1f}\n")

print(f"Detailed results saved to: {results_file}")
print()

# Side-by-side comparison
print("=" * 80)
print("SIDE-BY-SIDE COMPARISON")
print("=" * 80)
print()

# Load Infernal results for comparison
infernal_by_type = {}
with open(infernal_dir / 'diverse_test_results.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            rna_type = parts[1]
            success = parts[3] == 'True'
            if rna_type not in infernal_by_type:
                infernal_by_type[rna_type] = {'total': 0, 'success': 0}
            infernal_by_type[rna_type]['total'] += 1
            if success:
                infernal_by_type[rna_type]['success'] += 1

print(f"{'RNA Type':<15} {'Infernal+Rfam':<20} {'BLAST+RNAcentral':<20}")
print("-" * 60)

for rna_type in sorted(set(list(by_type.keys()) + list(infernal_by_type.keys()))):
    # BLAST results
    if rna_type in by_type:
        blast_total = len(by_type[rna_type])
        blast_success = sum(1 for r in by_type[rna_type] if r['success'])
        blast_rate = blast_success / blast_total * 100
        blast_str = f"{blast_success}/{blast_total} ({blast_rate:.0f}%)"
    else:
        blast_str = "N/A"
    
    # Infernal results
    if rna_type in infernal_by_type:
        inf_total = infernal_by_type[rna_type]['total']
        inf_success = infernal_by_type[rna_type]['success']
        inf_rate = inf_success / inf_total * 100
        inf_str = f"{inf_success}/{inf_total} ({inf_rate:.0f}%)"
    else:
        inf_str = "N/A"
    
    print(f"{rna_type:<15} {inf_str:<20} {blast_str:<20}")

print()
print("Overall:")
inf_total = sum(d['total'] for d in infernal_by_type.values())
inf_success = sum(d['success'] for d in infernal_by_type.values())
inf_rate = inf_success / inf_total * 100

blast_total = len(results)
blast_success = sum(1 for r in results if r['success'])
blast_rate = blast_success / blast_total * 100

print(f"{'TOTAL':<15} {inf_success}/{inf_total} ({inf_rate:.1f}%){'':<7} {blast_success}/{blast_total} ({blast_rate:.1f}%)")
