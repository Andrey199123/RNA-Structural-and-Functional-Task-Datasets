#!/usr/bin/env python3
"""
Generate balanced negative samples for RNA-protein interaction prediction.

Takes the positive pairs from RPI7317 and generates an equal number of
negative pairs by randomly pairing RNAs and proteins that don't interact.

Only uses RNAs and proteins that have structural features available.
"""

import csv
import random
from pathlib import Path

def main():
    # Get RNAs and proteins with features
    rna_features = set(
        p.stem.replace('_features', '') 
        for p in Path('RNA-Protein/data/structural_features/rna_features').glob('*.json')
    )
    prot_features = set(
        p.stem.replace('_features', '') 
        for p in Path('RNA-Protein/data/structural_features/protein_features').glob('*.json')
    )
    
    print(f"RNAs with features: {len(rna_features)}")
    print(f"Proteins with features: {len(prot_features)}")
    
    # Load positive pairs (only those with features)
    positive_pairs = set()
    with open('RNA-Protein/data/RPI7317/rpi7317_pairs.csv', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rna = row['RNA names']
            prot = row['Protein names']
            if rna in rna_features and prot in prot_features:
                positive_pairs.add((rna, prot))
    
    print(f"Positive pairs with features: {len(positive_pairs)}")
    
    # Generate negative pairs
    # These are random RNA-protein combinations that are NOT in the positive set
    rna_list = list(rna_features)
    prot_list = list(prot_features)
    
    negative_pairs = set()
    max_attempts = len(positive_pairs) * 100  # Prevent infinite loop
    attempts = 0
    
    while len(negative_pairs) < len(positive_pairs) and attempts < max_attempts:
        rna = random.choice(rna_list)
        prot = random.choice(prot_list)
        pair = (rna, prot)
        
        if pair not in positive_pairs and pair not in negative_pairs:
            negative_pairs.add(pair)
        
        attempts += 1
    
    print(f"Negative pairs generated: {len(negative_pairs)}")
    
    # Combine and save
    output_file = Path('RNA-Protein/data/RPI7317/rpi7317_pairs_balanced.csv')
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['RNA names', 'Protein names', 'Labels'])
        
        # Write positives
        for rna, prot in positive_pairs:
            writer.writerow([rna, prot, 1])
        
        # Write negatives
        for rna, prot in negative_pairs:
            writer.writerow([rna, prot, 0])
    
    print(f"\nSaved to: {output_file}")
    print(f"Total pairs: {len(positive_pairs) + len(negative_pairs)}")
    print(f"  - Positive: {len(positive_pairs)}")
    print(f"  - Negative: {len(negative_pairs)}")

if __name__ == "__main__":
    random.seed(42)  # For reproducibility
    main()
