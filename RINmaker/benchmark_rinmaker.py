#!/usr/bin/env python3
"""
Benchmark RINmaker performance on test proteins from the paper
Reproduces Table 2 results comparing RINmaker performance
"""

import os
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path

# Expected results from RINmaker paper (Table 2)
# CPU times in seconds (without hydrogen addition)
EXPECTED_PERFORMANCE = {
    '1al1': {'residues': 12, 'expected_time': 0.005},
    '1l2y': {'residues': 20, 'expected_time': 0.005},
    '1bk8': {'residues': 50, 'expected_time': 0.004},
    '1eod': {'residues': 100, 'expected_time': 0.009},
    '1hhq': {'residues': 150, 'expected_time': 0.009},
    '1ab5': {'residues': 250, 'expected_time': 0.015},
    '1bge': {'residues': 318, 'expected_time': 0.022},
    '1b1y': {'residues': 500, 'expected_time': 0.037},
    '2vji': {'residues': 598, 'expected_time': 0.049},
    '1miq': {'residues': 750, 'expected_time': 0.050},
    '1a3w': {'residues': 1000, 'expected_time': 0.061},
}

# Bond frequency from TOP2018 database (Table 1 in paper)
EXPECTED_BOND_FREQUENCY = {
    'VdW': 57.5,
    'H-Bond': 25.5,
    'Hydrophobic': 16.2,
    'Ionic': 0.5,
    'π-π stacking': 0.2,
    'π-cation': 0.1,
}

def run_rinmaker(pdb_file, rinmaker_binary='RINmaker/build/rinmaker', num_runs=5):
    """Run RINmaker on a PDB file and measure execution time"""
    output_file = pdb_file.replace('.pdb', '_rin.graphml')
    
    times = []
    
    for i in range(num_runs):
        start_time = time.time()
        
        try:
            result = subprocess.run(
                [rinmaker_binary, pdb_file, output_file],
                capture_output=True,
                text=True,
                timeout=60  # 1 minute timeout
            )
            
            elapsed_time = time.time() - start_time
            times.append(elapsed_time)
            
            if result.returncode != 0:
                print(f"✗ RINmaker failed: {result.stderr}")
                return None, None
                
        except subprocess.TimeoutExpired:
            print(f"✗ RINmaker timed out")
            return None, None
        except Exception as e:
            print(f"✗ Error running RINmaker: {e}")
            return None, None
    
    # Return average time and output file
    avg_time = sum(times) / len(times)
    return avg_time, output_file if os.path.exists(output_file) else None

def parse_graphml(graphml_file):
    """Parse GraphML output to extract bond statistics"""
    try:
        tree = ET.parse(graphml_file)
        root = tree.getroot()
        
        # Count edges (bonds) by type
        bond_counts = {
            'VdW': 0,
            'H-Bond': 0,
            'Hydrophobic': 0,
            'Ionic': 0,
            'π-π stacking': 0,
            'π-cation': 0,
        }
        
        # GraphML namespace
        ns = {'graphml': 'http://graphml.graphdrawing.org/xmlns'}
        
        # Count edges
        edges = root.findall('.//graphml:edge', ns)
        total_bonds = len(edges)
        
        # Parse edge attributes to count bond types
        for edge in edges:
            for data in edge.findall('graphml:data', ns):
                key = data.get('key')
                value = data.text
                
                # Map bond types (this depends on RINmaker's output format)
                if value:
                    if 'VdW' in value or 'vdw' in value.lower():
                        bond_counts['VdW'] += 1
                    elif 'HBond' in value or 'hbond' in value.lower():
                        bond_counts['H-Bond'] += 1
                    elif 'Hydrophobic' in value or 'hydrophobic' in value.lower():
                        bond_counts['Hydrophobic'] += 1
                    elif 'Ionic' in value or 'ionic' in value.lower():
                        bond_counts['Ionic'] += 1
                    elif 'PiPi' in value or 'pi-pi' in value.lower():
                        bond_counts['π-π stacking'] += 1
                    elif 'PiCation' in value or 'pi-cation' in value.lower():
                        bond_counts['π-cation'] += 1
        
        # Calculate frequencies
        frequencies = {}
        if total_bonds > 0:
            for bond_type, count in bond_counts.items():
                frequencies[bond_type] = (count / total_bonds) * 100
        
        return {
            'total_bonds': total_bonds,
            'bond_counts': bond_counts,
            'frequencies': frequencies
        }
        
    except Exception as e:
        print(f"✗ Error parsing GraphML: {e}")
        return None

def benchmark_structure(pdb_id, expected, rinmaker_binary):
    """Benchmark a single structure"""
    pdb_file = f'test_proteins/{pdb_id}.pdb'
    
    if not os.path.exists(pdb_file):
        print(f"✗ File not found: {pdb_file}")
        return None
    
    print(f"\nProcessing {pdb_id} ({expected['residues']} residues)...")
    print(f"Expected time: {expected['expected_time']:.3f}s")
    
    # Run RINmaker
    avg_time, output_file = run_rinmaker(pdb_file, rinmaker_binary)
    
    if avg_time is None:
        return None
    
    print(f"✓ Observed time: {avg_time:.3f}s")
    
    # Compare with expected
    speedup = expected['expected_time'] / avg_time if avg_time > 0 else 0
    if avg_time <= expected['expected_time'] * 2:  # Within 2x of expected
        print(f"✓ Performance: {speedup:.1f}x relative to paper")
    else:
        print(f"⚠ Performance: {speedup:.1f}x relative to paper (slower than expected)")
    
    # Parse output
    stats = None
    if output_file:
        stats = parse_graphml(output_file)
        if stats:
            print(f"  Total bonds: {stats['total_bonds']}")
    
    return {
        'pdb_id': pdb_id,
        'residues': expected['residues'],
        'expected_time': expected['expected_time'],
        'observed_time': avg_time,
        'speedup': speedup,
        'stats': stats
    }

def main():
    print("="*70)
    print("RINmaker Benchmark - Reproducing Paper Results")
    print("="*70)
    
    # Check if RINmaker binary exists
    rinmaker_binary = 'RINmaker/build/rinmaker'
    if not os.path.exists(rinmaker_binary):
        print(f"✗ RINmaker binary not found at {rinmaker_binary}")
        print("Please run install_rinmaker.sh first")
        return
    
    test_dir = 'test_proteins'
    if not os.path.exists(test_dir):
        print(f"✗ Test proteins directory not found: {test_dir}")
        print("Please run download_test_proteins.py first")
        return
    
    results = []
    
    for pdb_id, expected in EXPECTED_PERFORMANCE.items():
        result = benchmark_structure(pdb_id, expected, rinmaker_binary)
        if result:
            results.append(result)
    
    # Summary
    print("\n" + "="*70)
    print("BENCHMARK SUMMARY")
    print("="*70)
    print(f"{'PDB ID':<8} {'Residues':<10} {'Expected(s)':<12} {'Observed(s)':<12} {'Speedup':<10}")
    print("-"*70)
    
    for result in results:
        print(f"{result['pdb_id']:<8} {result['residues']:<10} "
              f"{result['expected_time']:<12.3f} {result['observed_time']:<12.3f} "
              f"{result['speedup']:<10.1f}x")
    
    # Calculate average speedup
    if results:
        avg_speedup = sum(r['speedup'] for r in results) / len(results)
        print(f"\nAverage speedup relative to paper: {avg_speedup:.1f}x")
        
        # Check if performance is acceptable
        if avg_speedup >= 0.5:  # Within 2x slower
            print("✓ Performance matches paper expectations")
        else:
            print("⚠ Performance is slower than expected")
    
    print("\n" + "="*70)
    print("Note: Performance may vary based on hardware and system load")
    print("Paper benchmarks were run on Intel i7-8565U CPU with 16GB RAM")
    print("="*70)

if __name__ == '__main__':
    main()
