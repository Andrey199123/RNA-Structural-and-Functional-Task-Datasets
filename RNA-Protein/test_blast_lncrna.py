import subprocess
from pathlib import Path
import time
import sys
import csv

# Check database
rnacentral_db = Path('data/rnacentral/rnacentral_active')

lncrnas = []
current_header = None
current_seq = []

with open('data/rnacentral/rnacentral_active.fasta', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_header and current_seq:
                seq = ''.join(current_seq)
                seq_len = len(seq)
                
                # Only LONG lncRNAs: 1000-2000nt
                if 1000 <= seq_len <= 2000:
                    header_lower = current_header.lower()
                    
                    if 'lncrna' in header_lower or 'long non-coding' in header_lower:
                        seq_id = current_header.split()[0][1:]
                        lncrnas.append({
                            'id': seq_id,
                            'seq': seq,
                            'length': seq_len,
                            'header': current_header
                        })
                        
                        if len(lncrnas) >= 30:
                            break
            
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)

print(f"Found {len(lncrnas)} long lncRNAs (>1000nt)")
print()

# Analyze length distribution
length_bins = {
    '1000-1200nt': 0,
    '1200-1400nt': 0,
    '1400-1600nt': 0,
    '1600-1800nt': 0,
    '1800-2000nt': 0
}

for rna in lncrnas:
    length = rna['length']
    if length < 1200:
        length_bins['1000-1200nt'] += 1
    elif length < 1400:
        length_bins['1200-1400nt'] += 1
    elif length < 1600:
        length_bins['1400-1600nt'] += 1
    elif length < 1800:
        length_bins['1600-1800nt'] += 1
    else:
        length_bins['1800-2000nt'] += 1

print("Length distribution:")
for bin_name, count in length_bins.items():
    print(f"  {bin_name}: {count} lncRNAs")
print()

# Create output directory
output_dir = Path('data/test_blast_long_lncrna')
output_dir.mkdir(exist_ok=True, parents=True)

# Test each lncRNA
results = []
success_count = 0

print("Testing each long lncRNA with BLAST + RNAcentral...")
print("Using optimized settings: max 50 homologs, faster MAFFT")
print()

for i, rna in enumerate(lncrnas, 1):
    print(f"[{i}/{len(lncrnas)}] {rna['id']} ({rna['length']} nt)")
    print(f"  Running BLAST...", end='', flush=True)
    
    # Write query file
    query_file = output_dir / f"{rna['id']}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna['id']}\n{rna['seq']}\n")
    
    try:
        start = time.time()
        
        # Step 1: BLAST search (limit to 50 for long RNAs to avoid timeout)
        blast_out = output_dir / f"{rna['id']}_blast.txt"
        cmd = [
            '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/blastn',
            '-query', str(query_file),
            '-db', 'data/rnacentral/rnacentral_active',
            '-evalue', '0.001',
            '-outfmt', '6 sseqid pident length',
            '-max_target_seqs', '50',  # Reduced from 100 for long RNAs
            '-num_threads', '4'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        
        # Parse BLAST results
        homolog_ids = set()
        identities = []
        
        with open(blast_out, 'w') as f:
            f.write(result.stdout)
        
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 3:
                    seq_id = parts[0]
                    pident = float(parts[1])
                    if pident >= 70:
                        homolog_ids.add(seq_id)
                        identities.append(pident)
        
        n_homologs = len(homolog_ids)
        avg_identity = sum(identities) / len(identities) if identities else 0
        
        print(f" {n_homologs} homologs", end='', flush=True)
        
        msa_generated = False
        a3m_file = None
        
        if n_homologs >= 10:
            # Step 2: Extract sequences
            print(f", extracting...", end='', flush=True)
            
            id_file = output_dir / f"{rna['id']}_ids.txt"
            with open(id_file, 'w') as f:
                for seq_id in homolog_ids:
                    f.write(f"{seq_id}\n")
            
            homologs_file = output_dir / f"{rna['id']}_homologs.fa"
            cmd = [
                '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/blastdbcmd',
                '-db', 'data/rnacentral/rnacentral_active',
                '-entry_batch', str(id_file),
                '-line_length', '100000',
                '-out', str(homologs_file)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            if homologs_file.exists() and homologs_file.stat().st_size > 0:
                # Step 3: Align with MAFFT (fast mode for long sequences)
                print(f" aligning...", end='', flush=True)
                
                aligned_file = output_dir / f"{rna['id']}_aligned.fa"
                # Use --retree 1 for faster alignment of long sequences
                cmd = [
                    '/Users/andreyvasilyev/miniconda3/envs/rna_msa/bin/mafft',
                    '--retree', '1',  # Faster mode
                    '--maxiterate', '0',  # Skip refinement
                    '--thread', '4',
                    str(homologs_file)
                ]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)  # 5 min timeout
                
                with open(aligned_file, 'w') as f:
                    f.write(result.stdout)
                
                if aligned_file.exists() and aligned_file.stat().st_size > 0:
                    # Step 4: Convert to A3M
                    print(f" converting to A3M...", end='', flush=True)
                    
                    sys.path.insert(0, 'src')
                    from rna_msa.alignment_formats import fasta_to_a3m, validate_a3m
                    
                    a3m_file = output_dir / f"{rna['id']}_aligned.a3m"
                    fasta_to_a3m(aligned_file, a3m_file, query_id=rna['id'])
                    
                    # Validate
                    stats = validate_a3m(a3m_file)
                    if stats['valid']:
                        msa_generated = True
                        success_count += 1
                        print(f" A3M generated ({stats['n_sequences']} seqs)")
                    else:
                        print(f" A3M invalid")
                else:
                    print(f" Alignment failed")
            else:
                print(f" Extraction failed")
        else:
            print(f" Insufficient homologs (need ≥10)")
        
        elapsed = time.time() - start
        print(f"  Time: {elapsed:.1f}s")
        
        results.append({
            'rna_id': rna['id'],
            'length': rna['length'],
            'msa_generated': msa_generated,
            'n_homologs': n_homologs,
            'avg_identity': avg_identity,
            'a3m_file': str(a3m_file) if a3m_file else None,
            'time': elapsed
        })
        
    except subprocess.TimeoutExpired:
        print(f" ✗ Timeout")
        results.append({
            'rna_id': rna['id'],
            'length': rna['length'],
            'msa_generated': False,
            'n_homologs': 0,
            'avg_identity': 0,
            'a3m_file': None,
            'time': 300
        })
    except Exception as e:
        print(f" Error: {str(e)[:50]}")
        results.append({
            'rna_id': rna['id'],
            'length': rna['length'],
            'msa_generated': False,
            'n_homologs': 0,
            'avg_identity': 0,
            'a3m_file': None,
            'time': 0
        })

# Save results

print(f"Total long lncRNAs tested: {len(lncrnas)}")
print(f"Successful MSAs: {success_count}")
print(f"Success rate: {success_count}/{len(lncrnas)} = {success_count/len(lncrnas)*100:.1f}%")


# By length bin
print("Success rate by length:")
for bin_name in ['1000-1200nt', '1200-1400nt', '1400-1600nt', '1600-1800nt', '1800-2000nt']:
    if bin_name == '1000-1200nt':
        bin_results = [r for r in results if r['length'] < 1200]
    elif bin_name == '1200-1400nt':
        bin_results = [r for r in results if 1200 <= r['length'] < 1400]
    elif bin_name == '1400-1600nt':
        bin_results = [r for r in results if 1400 <= r['length'] < 1600]
    elif bin_name == '1600-1800nt':
        bin_results = [r for r in results if 1600 <= r['length'] < 1800]
    else:
        bin_results = [r for r in results if r['length'] >= 1800]
    
    if bin_results:
        bin_success = sum(1 for r in bin_results if r['msa_generated'])
        print(f"  {bin_name}: {bin_success}/{len(bin_results)} ({bin_success/len(bin_results)*100:.1f}%)")


# Save to TSV
tsv_file = output_dir / 'long_lncrna_results.tsv'
with open(tsv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['rna_id', 'length', 'msa_generated', 'n_homologs', 'avg_identity', 'a3m_file', 'time'], delimiter='\t')
    writer.writeheader()
    writer.writerows(results)

print(f"Results saved to: {tsv_file}")

print("KEY FINDINGS:")
print(f"1. Long lncRNA (>1000nt) success rate: {success_count/len(lncrnas)*100:.1f}%")
print(f"2. Average homologs found: {sum(r['n_homologs'] for r in results)/len(results):.1f}")
