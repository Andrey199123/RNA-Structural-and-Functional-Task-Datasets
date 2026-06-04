from Bio import SeqIO
from collections import Counter, defaultdict
import json
from pathlib import Path
import re

def extract_rna_type(description):
    """Extract RNA type from RNAcentral header."""
    desc_lower = description.lower()
    
    # Check for specific RNA types in order of specificity
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

def extract_organism(description):
    """Extract organism from RNAcentral header if present."""
    # RNAcentral format often has "from X species" or organism name
    match = re.search(r'from (\d+) species', description)
    if match:
        return f"{match.group(1)} species"
    
    # Try to extract organism name patterns
    # This is approximate since RNAcentral headers vary
    if 'homo sapiens' in description.lower() or 'human' in description.lower():
        return 'Homo sapiens'
    elif 'mus musculus' in description.lower() or 'mouse' in description.lower():
        return 'Mus musculus'
    elif 'drosophila' in description.lower():
        return 'Drosophila'
    elif 'caenorhabditis elegans' in description.lower():
        return 'C. elegans'
    elif 'saccharomyces cerevisiae' in description.lower() or 'yeast' in description.lower():
        return 'S. cerevisiae'
    elif 'escherichia coli' in description.lower() or 'e. coli' in description.lower():
        return 'E. coli'
    elif 'arabidopsis' in description.lower():
        return 'Arabidopsis'
    else:
        return 'other/multiple'

def analyze_dataset(fasta_file, name):
    """Analyze RNA types and organisms in a dataset."""
    print(f"\nAnalyzing {name}...")
    
    rna_types = Counter()
    organisms = Counter()
    lengths = []
    total = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        total += 1
        rna_type = extract_rna_type(record.description)
        organism = extract_organism(record.description)
        
        rna_types[rna_type] += 1
        organisms[organism] += 1
        lengths.append(len(record.seq))
        
        if total % 100000 == 0:
            print(f"  Processed {total:,} sequences...")
    
    print(f"  Total: {total:,} sequences")
    
    return {
        'total': total,
        'rna_types': dict(rna_types),
        'organisms': dict(organisms),
        'length_stats': {
            'min': min(lengths),
            'max': max(lengths),
            'mean': sum(lengths) / len(lengths)
        }
    }

def generate_report(full_stats, subset_stats, output_file):
    """Generate comparison report."""
    report = []
    report.append("=" * 80)
    report.append("RNA CLASS DISTRIBUTION ANALYSIS")
    report.append("RNAcentral Full Dataset vs Subset")
    report.append("=" * 80)
    report.append("")
    
    # Overall counts
    report.append("OVERALL STATISTICS")
    report.append("-" * 80)
    report.append(f"Full dataset (filtered 50-2000 nt): {full_stats['total']:,} sequences")
    report.append(f"Subset (after sampling + dedup):    {subset_stats['total']:,} sequences")
    report.append(f"Subset percentage:                  {(subset_stats['total']/full_stats['total']*100):.2f}%")
    report.append("")
    
    # RNA type distribution
    report.append("RNA TYPE DISTRIBUTION")
    report.append("-" * 80)
    report.append(f"{'RNA Type':<20} {'Full Dataset':<20} {'Subset':<20} {'Representation':<15}")
    report.append("-" * 80)
    
    # Sort by full dataset count
    sorted_types = sorted(full_stats['rna_types'].items(), key=lambda x: x[1], reverse=True)
    
    for rna_type, full_count in sorted_types:
        subset_count = subset_stats['rna_types'].get(rna_type, 0)
        full_pct = (full_count / full_stats['total']) * 100
        subset_pct = (subset_count / subset_stats['total']) * 100 if subset_count > 0 else 0
        
        # Calculate representation ratio
        if full_pct > 0:
            representation = subset_pct / full_pct
        else:
            representation = 0
        
        report.append(f"{rna_type:<20} {full_count:>10,} ({full_pct:>5.2f}%)  {subset_count:>10,} ({subset_pct:>5.2f}%)  {representation:>6.2f}x")
    
    report.append("")
    
    # Organism distribution (top 10)
    report.append("ORGANISM DISTRIBUTION (Top 10)")
    report.append("-" * 80)
    report.append(f"{'Organism':<30} {'Full Dataset':<20} {'Subset':<20}")
    report.append("-" * 80)
    
    sorted_orgs = sorted(full_stats['organisms'].items(), key=lambda x: x[1], reverse=True)[:10]
    
    for organism, full_count in sorted_orgs:
        subset_count = subset_stats['organisms'].get(organism, 0)
        full_pct = (full_count / full_stats['total']) * 100
        subset_pct = (subset_count / subset_stats['total']) * 100 if subset_count > 0 else 0
        
        report.append(f"{organism:<30} {full_count:>10,} ({full_pct:>5.2f}%)  {subset_count:>10,} ({subset_pct:>5.2f}%)")
    
    report.append("")
    report.append("=" * 80)
    report.append("")
    report.append("INTERPRETATION")
    report.append("-" * 80)
    report.append("Representation ratio shows how well each RNA type is represented in the subset")
    report.append("compared to the full dataset:")
    report.append("  - 1.0x = perfectly proportional representation")
    report.append("  - >1.0x = overrepresented in subset")
    report.append("  - <1.0x = underrepresented in subset")
    report.append("")
    report.append("Random sampling should produce ratios close to 1.0x for all classes.")
    report.append("Significant deviations may indicate sampling bias or class imbalance issues.")
    report.append("=" * 80)
    
    report_text = "\n".join(report)
    
    with open(output_file, 'w') as f:
        f.write(report_text)
    
    print(f"\nReport saved to: {output_file}")
    print("\n" + report_text)
    
    return report_text

def main():
    # Paths
    full_dataset = Path("data/rnacentral/rnacentral_active.fasta")
    subset_dataset = Path("data/rnacentral_subset_balanced/final_subset.fasta")
    
    if not full_dataset.exists():
        print(f"Error: Full dataset not found at {full_dataset}")
        print("Please provide the path to the full RNAcentral FASTA file.")
        return
    
    if not subset_dataset.exists():
        print(f"Error: Subset not found at {subset_dataset}")
        return
    
    # Analyze both datasets
    print("Starting analysis...")
    print("This may take 20-30 minutes for the full RNAcentral dataset (40.7M sequences)...")
    
    full_stats = analyze_dataset(full_dataset, "Full RNAcentral")
    subset_stats = analyze_dataset(subset_dataset, "Balanced Subset")
    
    # Save JSON
    json_output = Path("data/rnacentral_subset_balanced/class_distribution.json")
    with open(json_output, 'w') as f:
        json.dump({
            'full_dataset': full_stats,
            'subset': subset_stats
        }, f, indent=2)
    print(f"\nJSON data saved to: {json_output}")
    
    # Generate report
    report_output = Path("data/rnacentral_subset_balanced/CLASS_DISTRIBUTION_REPORT.txt")
    generate_report(full_stats, subset_stats, report_output)

if __name__ == "__main__":
    main()
