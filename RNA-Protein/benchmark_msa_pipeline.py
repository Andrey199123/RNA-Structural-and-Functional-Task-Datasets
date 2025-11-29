#!/usr/bin/env python3
"""
Tests the MSA pipeline on diverse RNA sequences to see how well it works.
"""

import subprocess
from pathlib import Path
import sys
import time

print("Extracting diverse RNA sample from RNAcentral...")
print("Getting RNAs of different lengths and types...\n")

cmd = """
head -50000 data/rnacentral/rnacentral_active.fasta | \
awk '
/^>/ {
    if (seq && length(seq) >= 100 && length(seq) <= 2000) {
        print header;
        print seq;
        count++;
        if (count >= 50) exit;
    }
    header = $0;
    seq = "";
}
!/^>/ {
    seq = seq $0;
}
END {
    if (seq && length(seq) >= 100 && length(seq) <= 2000 && count < 50) {
        print header;
        print seq;
    }
}'
"""

result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
lines = result.stdout.strip().split('\n')

test_rnas = []
for i in range(0, len(lines), 2):
    if i+1 < len(lines):
        header = lines[i]
        seq = lines[i+1]
        seq_id = header.split()[0][1:]
        rna_type = "unknown"
        if "rRNA" in header:
            rna_type = "rRNA"
        elif "tRNA" in header:
            rna_type = "tRNA"
        elif "lncRNA" in header:
            rna_type = "lncRNA"
        elif "miRNA" in header:
            rna_type = "miRNA"
        elif "snoRNA" in header:
            rna_type = "snoRNA"
        
        test_rnas.append({
            'id': seq_id,
            'seq': seq,
            'length': len(seq),
            'type': rna_type,
            'header': header
        })

print(f"Selected {len(test_rnas)} test RNAs")
print(f"Length range: {min(r['length'] for r in test_rnas)}-{max(r['length'] for r in test_rnas)} nt")
print()

from collections import Counter
type_counts = Counter(r['type'] for r in test_rnas)
print("RNA types:")
for rna_type, count in type_counts.most_common():
    print(f"  {rna_type}: {count}")
print()

output_dir = Path('data/msa_benchmark')
output_dir.mkdir(exist_ok=True)

results = []
print("Testing MSA generation...")

for i, rna in enumerate(test_rnas[:30], 1):
    print(f"\n[{i}/30] {rna['id']} ({rna['type']}, {rna['length']} nt)")
    
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    start_time = time.time()
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
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60)
        
        hits = []
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                seq_id = parts[0]
                identity = float(parts[1])
                align_len = int(parts[2])
                
                if identity >= 70 and align_len >= rna['length'] * 0.5:
                    hits.append((seq_id, identity, align_len))
        
        search_time = time.time() - start_time
        
        n_homologs = len(hits)
        if n_homologs >= 20:
            quality = "Excellent"
        elif n_homologs >= 10:
            quality = "Good"
        elif n_homologs >= 5:
            quality = "Fair"
        elif n_homologs >= 1:
            quality = "Poor"
        else:
            quality = "None"
        
        avg_identity = sum(h[1] for h in hits) / len(hits) if hits else 0
        
        result_data = {
            'id': rna['id'],
            'type': rna['type'],
            'length': rna['length'],
            'n_homologs': n_homologs,
            'quality': quality,
            'avg_identity': avg_identity,
            'search_time': search_time,
            'success': n_homologs >= 5
        }
        results.append(result_data)
        
        print(f"  Homologs: {n_homologs} (avg identity: {avg_identity:.1f}%)")
        print(f"  Quality: {quality}")
        print(f"  Time: {search_time:.1f}s")
        
    except subprocess.TimeoutExpired:
        print(f"  Timeout")
        results.append({
            'id': rna['id'],
            'type': rna['type'],
            'length': rna['length'],
            'n_homologs': 0,
            'quality': 'Timeout',
            'avg_identity': 0,
            'search_time': 60,
            'success': False
        })
    except Exception as e:
        print(f"  Error: {e}")
        results.append({
            'id': rna['id'],
            'type': rna['type'],
            'length': rna['length'],
            'n_homologs': 0,
            'quality': 'Error',
            'avg_identity': 0,
            'search_time': 0,
            'success': False
        })

print("\nBENCHMARK RESULTS")

# Overall statistics
total = len(results)
successful = sum(1 for r in results if r['success'])
success_rate = successful / total * 100

print(f"\nOverall Success Rate: {successful}/{total} ({success_rate:.1f}%)")
print(f"  (Success = >=5 high-quality homologs)")

# By quality
quality_counts = Counter(r['quality'] for r in results)
print(f"\nMSA Quality Distribution:")
for quality in ['Excellent', 'Good', 'Fair', 'Poor', 'None', 'Timeout', 'Error']:
    count = quality_counts.get(quality, 0)
    pct = count / total * 100
    print(f"  {quality:12s}: {count:2d} ({pct:5.1f}%)")

# By RNA type
print(f"\nSuccess Rate by RNA Type:")
type_results = {}
for rna_type in set(r['type'] for r in results):
    type_data = [r for r in results if r['type'] == rna_type]
    type_success = sum(1 for r in type_data if r['success'])
    type_total = len(type_data)
    type_rate = type_success / type_total * 100 if type_total > 0 else 0
    type_results[rna_type] = (type_success, type_total, type_rate)
    print(f"  {rna_type:12s}: {type_success}/{type_total} ({type_rate:.1f}%)")

# MSA depth statistics
successful_results = [r for r in results if r['success']]
if successful_results:
    homolog_counts = [r['n_homologs'] for r in successful_results]
    identities = [r['avg_identity'] for r in successful_results]
    
    print(f"\nMSA Depth (for successful cases):")
    print(f"  Min:     {min(homolog_counts)} sequences")
    print(f"  Max:     {max(homolog_counts)} sequences")
    print(f"  Average: {sum(homolog_counts)/len(homolog_counts):.1f} sequences")
    print(f"  Median:  {sorted(homolog_counts)[len(homolog_counts)//2]} sequences")
    
    print(f"\nSequence Identity (for successful cases):")
    print(f"  Min:     {min(identities):.1f}%")
    print(f"  Max:     {max(identities):.1f}%")
    print(f"  Average: {sum(identities)/len(identities):.1f}%")

# Performance
avg_time = sum(r['search_time'] for r in results) / len(results)
print(f"\nPerformance:")
print(f"  Average search time: {avg_time:.1f}s per RNA")
print(f"  Total time: {sum(r['search_time'] for r in results):.1f}s")

# Length analysis
print(f"\nSuccess Rate by Length:")
short = [r for r in results if r['length'] < 300]
medium = [r for r in results if 300 <= r['length'] < 800]
long_rnas = [r for r in results if r['length'] >= 800]

for category, data, label in [
    (short, short, "Short (<300nt)"),
    (medium, medium, "Medium (300-800nt)"),
    (long_rnas, long_rnas, "Long (>=800nt)")
]:
    if data:
        cat_success = sum(1 for r in data if r['success'])
        cat_total = len(data)
        cat_rate = cat_success / cat_total * 100
        print(f"  {label:20s}: {cat_success}/{cat_total} ({cat_rate:.1f}%)")

print("\n" + "=" * 80)
print("\nCONCLUSION:")
print("=" * 80)

if success_rate >= 70:
    print("EXCELLENT: Pipeline works well for most RNAs")
    print("  Can reliably generate MSAs for RNA-protein interaction prediction")
elif success_rate >= 50:
    print("GOOD: Pipeline works for majority of RNAs")
    print("  Suitable for MSA-based feature extraction")
elif success_rate >= 30:
    print("FAIR: Pipeline works for some RNAs")
    print("  May need to filter dataset or use hybrid approach")
else:
    print("POOR: Pipeline struggles with most RNAs")
    print("  Need alternative strategies or different dataset")

print(f"\nKey findings:")
print(f"  {success_rate:.0f}% of RNAs can generate good MSAs (>=5 homologs)")
if successful_results:
    print(f"  Average MSA depth: {sum(r['n_homologs'] for r in successful_results)/len(successful_results):.0f} sequences")
    print(f"  Average identity: {sum(r['avg_identity'] for r in successful_results)/len(successful_results):.0f}%")
print(f"  Fast search: ~{avg_time:.0f}s per RNA")

# Save detailed results
results_file = output_dir / 'benchmark_results.tsv'
with open(results_file, 'w') as f:
    f.write("RNA_ID\tType\tLength\tHomologs\tQuality\tAvg_Identity\tSearch_Time\tSuccess\n")
    for r in results:
        f.write(f"{r['id']}\t{r['type']}\t{r['length']}\t{r['n_homologs']}\t{r['quality']}\t{r['avg_identity']:.1f}\t{r['search_time']:.1f}\t{r['success']}\n")

print(f"\nDetailed results saved to: {results_file}")
