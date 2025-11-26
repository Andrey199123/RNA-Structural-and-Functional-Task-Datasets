#!/usr/bin/env python3
"""
Comprehensive test of all RNA MSA strategies on real RPI7317 data
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from rna_msa.pipeline import RNASequence, RNAMSAPipeline, RNAStrategy
import json

# Initialize pipeline
pipeline = RNAMSAPipeline(
    rfam_db_path=Path("data/rfam/Rfam.cm"),
    rfam_seed_path=Path("data/rfam/Rfam.seed"),
    sequence_db_path=Path("data/rfam/rfam_sequences.fasta"),
    output_dir=Path("data/rna_msa_comprehensive_test"),
    n_cpu=2
)

# Load first 10 RNAs from RPI7317
sequences = RNASequence.from_fasta(Path("data/RPI7317/RPI7317_IncRNA_sequence.fasta"))[:10]

print("=" * 80)
print("COMPREHENSIVE RNA MSA STRATEGY TEST")
print("=" * 80)
print(f"\nTesting {len(sequences)} RNAs from RPI7317 dataset")
print(f"Strategies: Infernal+Rfam, RNAcmap-style, Domain Split, Single-Sequence\n")

results = {
    'infernal_rfam': {'tested': [], 'succeeded': [], 'failed': []},
    'rnacmap_style': {'tested': [], 'succeeded': [], 'failed': []},
    'domain_split': {'tested': [], 'succeeded': [], 'failed': []},
    'single_sequence': {'tested': [], 'succeeded': [], 'failed': []}
}

for rna in sequences:
    print(f"\n{'='*80}")
    print(f"RNA: {rna.id} ({rna.length} nt)")
    print(f"{'='*80}")
    
    # Determine which strategies to test based on length
    if rna.length > 1000:
        strategies_to_test = [RNAStrategy.DOMAIN_SPLIT]
        print(f"Long RNA (>{rna.length}nt) - Testing: Domain Split only")
    else:
        strategies_to_test = [RNAStrategy.INFERNAL_RFAM, RNAStrategy.RNACMAP_STYLE, RNAStrategy.SINGLE_SEQUENCE]
        print(f"Short RNA ({rna.length}nt) - Testing: Rfam, RNAcmap, Single-Sequence")
    
    for strategy in strategies_to_test:
        strategy_name = strategy.value
        results[strategy_name]['tested'].append(rna.id)
        
        print(f"\n  Testing {strategy_name}...")
        try:
            result = pipeline.process_rna(rna, strategy=strategy, force_strategy=True)
            
            if result['success']:
                results[strategy_name]['succeeded'].append(rna.id)
                print(f"    ✓ SUCCESS")
                if 'msa_depth' in result:
                    print(f"      MSA depth: {result['msa_depth']}")
                if 'rfam_family' in result:
                    print(f"      Rfam family: {result['rfam_family']}")
                if 'n_domains' in result:
                    print(f"      Domains: {result['n_domains']}")
            else:
                results[strategy_name]['failed'].append(rna.id)
                print(f"    ✗ FAILED: {result.get('error', 'Unknown')}")
                
        except Exception as e:
            results[strategy_name]['failed'].append(rna.id)
            print(f"    ✗ EXCEPTION: {str(e)[:100]}")

# Summary
print(f"\n\n{'='*80}")
print("SUMMARY")
print(f"{'='*80}\n")

for strategy_name, data in results.items():
    if not data['tested']:
        continue
    
    total = len(data['tested'])
    succeeded = len(data['succeeded'])
    failed = len(data['failed'])
    success_rate = (succeeded / total * 100) if total > 0 else 0
    
    print(f"{strategy_name.upper()}")
    print(f"  Tested: {total}")
    print(f"  Succeeded: {succeeded}")
    print(f"  Failed: {failed}")
    print(f"  Success rate: {success_rate:.1f}%")
    if data['succeeded']:
        print(f"  Successful RNAs: {', '.join(data['succeeded'][:5])}")
    if data['failed']:
        print(f"  Failed RNAs: {', '.join(data['failed'][:5])}")
    print()

# Save results
with open('data/rna_msa_comprehensive_test/test_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"Detailed results saved to: data/rna_msa_comprehensive_test/test_results.json")
