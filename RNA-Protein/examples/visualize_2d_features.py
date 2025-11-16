#!/usr/bin/env python3
"""
Visualize 2D RNA structural features

This script creates visualizations of the extracted features to help
understand what the pipeline is generating.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse


def plot_rna_features(rna_id, features, output_dir=None):
    """
    Create comprehensive visualization of RNA features
    
    Args:
        rna_id: RNA identifier
        features: Dictionary with conservation, coevolution, accessibility
        output_dir: Optional directory to save plots
    """
    conservation = features['conservation']
    coevolution = features['coevolution']
    accessibility = features['accessibility']
    msa_depth = features.get('msa_depth', 0)
    
    length = len(conservation)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # Title
    fig.suptitle(f'RNA: {rna_id} (Length: {length} nt, MSA depth: {msa_depth})', 
                 fontsize=14, fontweight='bold')
    
    # 1. Conservation profile
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(conservation, linewidth=1.5, color='#2E86AB')
    ax1.fill_between(range(length), conservation, alpha=0.3, color='#2E86AB')
    ax1.set_xlabel('Position')
    ax1.set_ylabel('Conservation Score')
    ax1.set_title('Conservation Profile (higher = more conserved)')
    ax1.set_ylim([0, 1])
    ax1.grid(True, alpha=0.3)
    
    # 2. Accessibility profile
    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(accessibility, linewidth=1.5, color='#A23B72')
    ax2.fill_between(range(length), accessibility, alpha=0.3, color='#A23B72')
    ax2.set_xlabel('Position')
    ax2.set_ylabel('Accessibility Score')
    ax2.set_title('Accessibility Profile (higher = more accessible/unpaired)')
    ax2.set_ylim([0, 1])
    ax2.grid(True, alpha=0.3)
    
    # 3. Coevolution/Contact map
    ax3 = fig.add_subplot(gs[2, 0])
    im = ax3.imshow(coevolution, cmap='YlOrRd', aspect='auto', interpolation='nearest')
    ax3.set_xlabel('Position')
    ax3.set_ylabel('Position')
    ax3.set_title('Coevolution/Contact Map')
    plt.colorbar(im, ax=ax3, label='Coevolution Score')
    
    # 4. Combined feature summary
    ax4 = fig.add_subplot(gs[2, 1])
    
    # Identify interesting regions
    high_conservation = conservation > 0.7
    high_accessibility = accessibility > 0.7
    potential_binding = high_accessibility & (conservation > 0.5)
    
    # Plot regions
    x = np.arange(length)
    ax4.fill_between(x, 0, 1, where=high_conservation, alpha=0.3, 
                     color='blue', label='Highly conserved')
    ax4.fill_between(x, 0, 1, where=high_accessibility, alpha=0.3, 
                     color='red', label='Highly accessible')
    ax4.fill_between(x, 0, 1, where=potential_binding, alpha=0.5, 
                     color='green', label='Potential binding site')
    
    ax4.set_xlabel('Position')
    ax4.set_ylabel('Feature Presence')
    ax4.set_title('Feature Summary')
    ax4.set_ylim([0, 1])
    ax4.legend(loc='upper right')
    ax4.grid(True, alpha=0.3)
    
    # Save or show
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{rna_id}_features.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {output_file}")
    else:
        plt.show()
    
    plt.close()


def plot_feature_statistics(features_file, output_dir=None):
    """
    Plot statistics across all RNAs in the dataset
    
    Args:
        features_file: NPZ file with features
        output_dir: Optional directory to save plots
    """
    data = np.load(features_file, allow_pickle=True)
    
    # Collect statistics
    lengths = []
    msa_depths = []
    avg_conservation = []
    avg_accessibility = []
    avg_coevolution = []
    
    for rna_id in data.files:
        features = data[rna_id].item()
        
        conservation = features['conservation']
        coevolution = features['coevolution']
        accessibility = features['accessibility']
        msa_depth = features.get('msa_depth', 0)
        
        lengths.append(len(conservation))
        msa_depths.append(msa_depth)
        avg_conservation.append(conservation.mean())
        avg_accessibility.append(accessibility.mean())
        avg_coevolution.append(coevolution.mean())
    
    # Create plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Feature Statistics Across {len(data.files)} RNAs', 
                 fontsize=14, fontweight='bold')
    
    # 1. RNA length distribution
    axes[0, 0].hist(lengths, bins=50, color='#2E86AB', alpha=0.7, edgecolor='black')
    axes[0, 0].set_xlabel('RNA Length (nt)')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('RNA Length Distribution')
    axes[0, 0].axvline(np.median(lengths), color='red', linestyle='--', 
                       label=f'Median: {np.median(lengths):.0f}')
    axes[0, 0].legend()
    
    # 2. MSA depth distribution
    axes[0, 1].hist(msa_depths, bins=50, color='#F18F01', alpha=0.7, edgecolor='black')
    axes[0, 1].set_xlabel('MSA Depth')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('MSA Depth Distribution')
    axes[0, 1].axvline(np.median(msa_depths), color='red', linestyle='--', 
                       label=f'Median: {np.median(msa_depths):.0f}')
    axes[0, 1].legend()
    
    # 3. Conservation distribution
    axes[0, 2].hist(avg_conservation, bins=50, color='#2E86AB', alpha=0.7, edgecolor='black')
    axes[0, 2].set_xlabel('Average Conservation')
    axes[0, 2].set_ylabel('Count')
    axes[0, 2].set_title('Conservation Distribution')
    axes[0, 2].axvline(np.median(avg_conservation), color='red', linestyle='--', 
                       label=f'Median: {np.median(avg_conservation):.2f}')
    axes[0, 2].legend()
    
    # 4. Accessibility distribution
    axes[1, 0].hist(avg_accessibility, bins=50, color='#A23B72', alpha=0.7, edgecolor='black')
    axes[1, 0].set_xlabel('Average Accessibility')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('Accessibility Distribution')
    axes[1, 0].axvline(np.median(avg_accessibility), color='red', linestyle='--', 
                       label=f'Median: {np.median(avg_accessibility):.2f}')
    axes[1, 0].legend()
    
    # 5. Coevolution distribution
    axes[1, 1].hist(avg_coevolution, bins=50, color='#C73E1D', alpha=0.7, edgecolor='black')
    axes[1, 1].set_xlabel('Average Coevolution')
    axes[1, 1].set_ylabel('Count')
    axes[1, 1].set_title('Coevolution Distribution')
    axes[1, 1].axvline(np.median(avg_coevolution), color='red', linestyle='--', 
                       label=f'Median: {np.median(avg_coevolution):.2f}')
    axes[1, 1].legend()
    
    # 6. MSA depth vs RNA length
    axes[1, 2].scatter(lengths, msa_depths, alpha=0.5, s=20, color='#2E86AB')
    axes[1, 2].set_xlabel('RNA Length (nt)')
    axes[1, 2].set_ylabel('MSA Depth')
    axes[1, 2].set_title('MSA Depth vs RNA Length')
    axes[1, 2].set_xscale('log')
    axes[1, 2].set_yscale('log')
    
    plt.tight_layout()
    
    # Save or show
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / "feature_statistics.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved statistics plot to {output_file}")
    else:
        plt.show()
    
    plt.close()
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("Feature Statistics Summary")
    print("=" * 60)
    print(f"Number of RNAs: {len(data.files)}")
    print(f"\nRNA Length:")
    print(f"  Min: {min(lengths)} nt")
    print(f"  Max: {max(lengths)} nt")
    print(f"  Mean: {np.mean(lengths):.1f} nt")
    print(f"  Median: {np.median(lengths):.1f} nt")
    print(f"\nMSA Depth:")
    print(f"  Min: {min(msa_depths)}")
    print(f"  Max: {max(msa_depths)}")
    print(f"  Mean: {np.mean(msa_depths):.1f}")
    print(f"  Median: {np.median(msa_depths):.1f}")
    print(f"  Single-sequence (depth=0): {sum(1 for d in msa_depths if d == 0)}")
    print(f"\nAverage Conservation:")
    print(f"  Min: {min(avg_conservation):.3f}")
    print(f"  Max: {max(avg_conservation):.3f}")
    print(f"  Mean: {np.mean(avg_conservation):.3f}")
    print(f"\nAverage Accessibility:")
    print(f"  Min: {min(avg_accessibility):.3f}")
    print(f"  Max: {max(avg_accessibility):.3f}")
    print(f"  Mean: {np.mean(avg_accessibility):.3f}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Visualize 2D RNA structural features"
    )
    
    parser.add_argument(
        "--features_file",
        type=Path,
        default=Path("data/embeddings/rna_2d_features_rpi7317.npz"),
        help="NPZ file with 2D features"
    )
    
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("data/rna_msa_output/visualizations"),
        help="Output directory for plots"
    )
    
    parser.add_argument(
        "--rna_ids",
        nargs='+',
        default=None,
        help="Specific RNA IDs to visualize (default: first 5)"
    )
    
    parser.add_argument(
        "--n_examples",
        type=int,
        default=5,
        help="Number of example RNAs to visualize"
    )
    
    parser.add_argument(
        "--statistics_only",
        action="store_true",
        help="Only plot statistics, not individual RNAs"
    )
    
    args = parser.parse_args()
    
    # Check if features file exists
    if not args.features_file.exists():
        print(f"Features file not found: {args.features_file}")
        print("Please run the RNA MSA pipeline first:")
        print("  python run_rna_msa_pipeline.py --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta")
        return 1
    
    # Load features
    data = np.load(args.features_file, allow_pickle=True)
    print(f"Loaded features for {len(data.files)} RNAs")
    
    # Plot overall statistics
    print("\nGenerating statistics plots...")
    plot_feature_statistics(args.features_file, args.output_dir)
    
    if not args.statistics_only:
        # Select RNAs to visualize
        if args.rna_ids:
            rna_ids = args.rna_ids
        else:
            rna_ids = data.files[:args.n_examples]
        
        print(f"\nGenerating individual RNA plots for {len(rna_ids)} RNAs...")
        
        # Plot individual RNAs
        for i, rna_id in enumerate(rna_ids, 1):
            if rna_id not in data.files:
                print(f"Warning: RNA {rna_id} not found in features file")
                continue
            
            print(f"[{i}/{len(rna_ids)}] Plotting {rna_id}...")
            features = data[rna_id].item()
            plot_rna_features(rna_id, features, args.output_dir)
    
    print(f"\nAll plots saved to: {args.output_dir}")
    print("Done!")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
