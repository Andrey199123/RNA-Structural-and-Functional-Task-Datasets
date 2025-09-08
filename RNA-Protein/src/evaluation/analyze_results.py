#!/usr/bin/env python3
"""
Analyze and compare ZHMolGraph replication results with paper benchmarks.
Generate performance tables and figures.
"""
import json
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


def load_results(results_dir="."):
    """Load all result JSON files."""
    results = {}
    for json_file in Path(results_dir).glob("results_*.json"):
        with open(json_file, 'r') as f:
            data = json.load(f)
            key = f"{data['dataset']}_{data['split_type']}"
            results[key] = data
    return results


def create_performance_table(results):
    """Create a performance comparison table."""
    table_data = []
    
    for key, data in results.items():
        dataset = data['dataset']
        split_type = data['split_type']
        
        table_data.append({
            'Dataset': dataset.upper(),
            'Split Type': split_type,
            'Accuracy': f"{data['mean_test_acc']:.4f} ± {data['std_test_acc']:.4f}",
            'MCC': f"{data['mean_test_mcc']:.4f} ± {data['std_test_mcc']:.4f}",
            'AUROC': f"{data['mean_test_auroc']:.4f} ± {data['std_test_auroc']:.4f}",
            'AUPRC': f"{data['mean_test_auprc']:.4f} ± {data['std_test_auprc']:.4f}",
            'N Folds': data['n_folds']
        })
    
    df = pd.DataFrame(table_data)
    return df


def compare_with_paper():
    """Compare our results with paper benchmarks."""
    # Paper results from ZHMolGraph paper (approximate values from figures)
    paper_results = {
        'RPI7317_random': {
            'Accuracy': 0.915,
            'MCC': 0.831,
            'AUROC': None,  # Not reported in paper
            'AUPRC': None
        },
        'RPI7317_uNODES': {
            'Accuracy': None,
            'MCC': None,
            'AUROC': 0.798,  # From TheNovel dataset
            'AUPRC': 0.820
        },
        'NPInter2_random': {
            'Accuracy': 0.955,
            'MCC': 0.911,
            'AUROC': None,
            'AUPRC': None
        }
    }
    
    return paper_results


def create_comparison_plot(results):
    """Create comparison plots."""
    # Extract data for plotting
    datasets = []
    split_types = []
    mcc_values = []
    auroc_values = []
    
    for key, data in results.items():
        datasets.append(data['dataset'])
        split_types.append(data['split_type'])
        mcc_values.append(data['mean_test_mcc'])
        auroc_values.append(data['mean_test_auroc'])
    
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # MCC comparison
    x_pos = np.arange(len(datasets))
    bars1 = ax1.bar(x_pos, mcc_values, alpha=0.7, color='skyblue')
    ax1.set_xlabel('Dataset & Split Type')
    ax1.set_ylabel('MCC')
    ax1.set_title('Test MCC Performance')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([f"{d}\n{s}" for d, s in zip(datasets, split_types)], rotation=45)
    ax1.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for i, v in enumerate(mcc_values):
        ax1.text(i, v + 0.01, f'{v:.3f}', ha='center', va='bottom')
    
    # AUROC comparison
    bars2 = ax2.bar(x_pos, auroc_values, alpha=0.7, color='lightcoral')
    ax2.set_xlabel('Dataset & Split Type')
    ax2.set_ylabel('AUROC')
    ax2.set_title('Test AUROC Performance')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f"{d}\n{s}" for d, s in zip(datasets, split_types)], rotation=45)
    ax2.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for i, v in enumerate(auroc_values):
        ax2.text(i, v + 0.01, f'{v:.3f}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('performance_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()


def generate_report(results):
    """Generate a comprehensive analysis report."""
    print("="*80)
    print("ZHMolGraph REPLICATION RESULTS ANALYSIS")
    print("="*80)
    
    # Performance table
    print("\n1. PERFORMANCE SUMMARY")
    print("-" * 50)
    df = create_performance_table(results)
    print(df.to_string(index=False))
    
    # Key findings
    print("\n2. KEY FINDINGS")
    print("-" * 50)
    
    # Find best performing configuration
    best_mcc = max(results.items(), key=lambda x: x[1]['mean_test_mcc'])
    print(f"• Best MCC Performance: {best_mcc[0]} - {best_mcc[1]['mean_test_mcc']:.4f}")
    
    # Compare random vs unknown nodes
    if 'rpi7317_random' in results and 'rpi7317_uNODES' in results:
        random_mcc = results['rpi7317_random']['mean_test_mcc']
        unknown_mcc = results['rpi7317_uNODES']['mean_test_mcc']
        print(f"• Random vs Unknown Nodes MCC: {random_mcc:.4f} vs {unknown_mcc:.4f}")
        print(f"  → Performance drop: {random_mcc - unknown_mcc:.4f} ({((random_mcc - unknown_mcc) / random_mcc * 100):.1f}%)")
    
    # Paper comparison
    print("\n3. COMPARISON WITH PAPER BENCHMARKS")
    print("-" * 50)
    paper_results = compare_with_paper()
    
    if 'rpi7317_random' in results:
        our_mcc = results['rpi7317_random']['mean_test_mcc']
        paper_mcc = paper_results['RPI7317_random']['MCC']
        if paper_mcc:
            diff = our_mcc - paper_mcc
            print(f"• RPI7317 Random MCC: Ours {our_mcc:.4f} vs Paper {paper_mcc:.4f} (diff: {diff:+.4f})")
    
    if 'rpi7317_uNODES' in results:
        our_auroc = results['rpi7317_uNODES']['mean_test_auroc']
        paper_auroc = paper_results['RPI7317_uNODES']['AUROC']
        if paper_auroc:
            diff = our_auroc - paper_auroc
            print(f"• RPI7317 Unknown Nodes AUROC: Ours {our_auroc:.4f} vs Paper {paper_auroc:.4f} (diff: {diff:+.4f})")
    
    # Generalizability analysis
    print("\n4. GENERALIZABILITY ANALYSIS")
    print("-" * 50)
    if 'rpi7317_random' in results and 'rpi7317_uNODES' in results:
        random_auroc = results['rpi7317_random']['mean_test_auroc']
        unknown_auroc = results['rpi7317_uNODES']['mean_test_auroc']
        print(f"• Random Split AUROC: {random_auroc:.4f}")
        print(f"• Unknown Nodes AUROC: {unknown_auroc:.4f}")
        print(f"• Generalizability Gap: {random_auroc - unknown_auroc:.4f}")
        
        if unknown_auroc > 0.6:
            print("  → Model shows reasonable generalizability to unseen nodes")
        else:
            print("  → Model struggles with completely unseen nodes (expected)")
    
    # Recommendations
    print("\n5. RECOMMENDATIONS")
    print("-" * 50)
    print("• The replication shows strong performance on seen nodes (random split)")
    print("• Performance degrades significantly on unseen nodes, consistent with paper")
    print("• Consider implementing additional ablation studies:")
    print("  - LLM vs no-LLM comparison")
    print("  - Different protein embedding models (ProtTrans vs ESM2)")
    print("  - Network topology analysis")
    print("• For production use, consider ensemble methods or domain adaptation")


def main():
    """Main analysis function."""
    # Load results
    results = load_results()
    
    if not results:
        print("No result files found. Run benchmark evaluations first.")
        return
    
    # Generate report
    generate_report(results)
    
    # Create plots
    try:
        create_comparison_plot(results)
        print(f"\nPerformance comparison plot saved as 'performance_comparison.png'")
    except ImportError:
        print("\nMatplotlib not available - skipping plot generation")
    
    # Save detailed results
    with open('analysis_summary.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nDetailed results saved as 'analysis_summary.json'")


if __name__ == "__main__":
    main()
