import sys
from pathlib import Path
from Bio import SeqIO
import subprocess
import numpy as np
import argparse
from tqdm import tqdm
from collections import defaultdict
import re

cdhit_identity = 0.95
min_class_representation = 1000
dedup_factor = 1.6


def extract_rna_type(description):
    desc_lower = description.lower()
    
    if 'rrna' in desc_lower or 'ribosomal rna' in desc_lower or 'ribosomal' in desc_lower:
        return 'rRNA'
    elif 'trna' in desc_lower or 'transfer rna' in desc_lower:
        return 'tRNA'
    elif 'mirna' in desc_lower or 'microrna' in desc_lower or 'mir-' in desc_lower:
        return 'miRNA'
    elif 'snrna' in desc_lower or 'small nuclear' in desc_lower:
        return 'snRNA'
    elif 'snorna' in desc_lower or 'small nucleolar' in desc_lower:
        return 'snoRNA'
    elif 'lncrna' in desc_lower or 'long non-coding' in desc_lower or 'long noncoding' in desc_lower:
        return 'lncRNA'
    elif 'mrna' in desc_lower or 'messenger rna' in desc_lower:
        return 'mRNA'
    elif 'ncrna' in desc_lower or 'non-coding rna' in desc_lower or 'noncoding rna' in desc_lower:
        return 'ncRNA'
    elif 'riboswitch' in desc_lower:
        return 'riboswitch'
    elif 'ribozyme' in desc_lower:
        return 'ribozyme'
    elif 'telomerase' in desc_lower:
        return 'telomerase RNA'
    elif 'vault' in desc_lower:
        return 'vault RNA'
    elif 'y rna' in desc_lower or 'y_rna' in desc_lower:
        return 'Y RNA'
    elif 'srp' in desc_lower or 'signal recognition' in desc_lower:
        return 'SRP RNA'
    elif 'rnase' in desc_lower:
        return 'RNase RNA'
    elif 'antisense' in desc_lower:
        return 'antisense RNA'
    elif 'crispr' in desc_lower:
        return 'CRISPR RNA'
    else:
        return 'other'


def count_classes(input_file):
    class_counts = defaultdict(int)
    total = 0
    
    for record in tqdm(SeqIO.parse(input_file, "fasta"), desc="Counting sequences"):
        rna_type = extract_rna_type(record.description)
        class_counts[rna_type] += 1
        total += 1
    
    return class_counts, total


def calculate_balanced_targets(class_counts, total_count, target_total, min_per_class=min_class_representation):
    sample_total = int(target_total * dedup_factor)
    num_classes = len(class_counts)
    min_budget = min_per_class * num_classes
    
    if min_budget > sample_total:
        min_per_class = sample_total // num_classes
        min_budget = min_per_class * num_classes
    
    remaining_budget = sample_total - min_budget
    
    targets = {}
    for rna_type, count in class_counts.items():
        target = min_per_class
        proportion = count / total_count
        proportional_share = int(remaining_budget * proportion)
        target += proportional_share
        target = min(target, count)
        targets[rna_type] = target
    
    return targets


def balanced_sample(input_file, output_file, targets, class_counts):
    np.random.seed(42)
    class_indices = {}
    
    for rna_type, target in targets.items():
        available = class_counts[rna_type]
        if target > 0 and available > 0:
            actual_target = min(target, available)
            indices = set(np.random.choice(available, actual_target, replace=False))
            class_indices[rna_type] = indices
    
    class_current_idx = defaultdict(int)
    sampled = []
    
    for record in tqdm(SeqIO.parse(input_file, "fasta"), desc="Sampling sequences"):
        rna_type = extract_rna_type(record.description)
        current_idx = class_current_idx[rna_type]
        
        if rna_type in class_indices and current_idx in class_indices[rna_type]:
            sampled.append(record)
        
        class_current_idx[rna_type] += 1
    
    SeqIO.write(sampled, output_file, "fasta")
    return output_file


def deduplicate_sequences(input_file, output_file, num_threads=12):
    import time
    
    cdhit_path = 'cd-hit-est'
    try:
        result = subprocess.run(['which', 'cd-hit-est'], capture_output=True, text=True)
        if result.returncode == 0:
            cdhit_path = result.stdout.strip()
    except:
        pass
    
    input_count = sum(1 for _ in SeqIO.parse(input_file, "fasta"))
    
    cmd = [
        cdhit_path,
        '-i', str(input_file),
        '-o', str(output_file),
        '-c', str(cdhit_identity),
        '-n', '10',
        '-M', '0',
        '-T', str(num_threads),
        '-d', '0',
        '-g', '1',
        '-aS', '0.8',
        '-aL', '0.8'
    ]
    
    print(f"Running CD-HIT-EST on {input_count:,} sequences...")
    
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        
        start_time = time.time()
        last_progress = 0
        
        for line in process.stdout:
            match = re.search(r'(\d+)\s+finished', line)
            if match:
                current = int(match.group(1))
                if current > last_progress and current > 1000:
                    progress = (current / input_count) * 100
                    elapsed = time.time() - start_time
                    rate = current / elapsed
                    remaining = input_count - current
                    eta_seconds = remaining / rate if rate > 0 else 0
                    eta_mins = int(eta_seconds // 60)
                    print(f"Progress: {progress:.1f}% | ETA: {eta_mins}m", end='\r')
                    last_progress = current
        
        process.wait()
        
        if process.returncode != 0:
            print(f"\nError: CD-HIT-EST failed")
            sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)
    
    output_count = sum(1 for _ in SeqIO.parse(output_file, "fasta"))
    print(f"\nDone: {output_count:,} sequences after deduplication")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(description="Create balanced RNAcentral subset")
    parser.add_argument('--input', type=str, default='data/rnacentral/rnacentral_active.fasta')
    parser.add_argument('--output-dir', type=str, default='data/rnacentral_subset_balanced')
    parser.add_argument('--target-size', type=int, default=400000)
    parser.add_argument('--min-per-class', type=int, default=1000)
    parser.add_argument('--threads', type=int, default=12)
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)
    
    print("Balanced RNAcentral Subset Creation")
    print(f"Target: {args.target_size:,} sequences\n")
    
    class_counts, total_count = count_classes(input_path)
    targets = calculate_balanced_targets(class_counts, total_count, args.target_size, args.min_per_class)
    
    sampled_file = output_dir / "01_balanced_sampled.fasta"
    balanced_sample(input_path, sampled_file, targets, class_counts)
    
    final_file = output_dir / "final_subset.fasta"
    deduplicate_sequences(sampled_file, final_file, args.threads)
    
    print(f"\nComplete: {final_file}")


if __name__ == "__main__":
    main()
