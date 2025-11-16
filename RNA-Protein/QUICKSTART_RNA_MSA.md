# Quick Start Guide: RNA MSA Pipeline

## Step 1: Install Dependencies

```bash
# Install Python packages
pip install -r requirements_rna_msa.txt

# Install external tools (Ubuntu/Debian)
sudo apt-get install infernal vienna-rna mafft

# Or on macOS
brew install infernal viennarna mafft
```

## Step 2: Check Installation

```bash
python run_rna_msa_pipeline.py --check_deps
```

Expected output:
```
Checking dependencies...
  ✓ infernal
  ✓ mafft
  ✓ RNAfold
  ✓ python_RNA
```

## Step 3: Test Pipeline (Small Scale)

Run on first 5 sequences to verify everything works:

```bash
python run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --output_dir data/rna_msa_output \
    --test_mode \
    --verbose
```

This should complete in a few minutes and create:
- `data/rna_msa_output/` directory with MSA files
- `data/rna_msa_output/pipeline_results.json` with summary
- Log file `rna_msa_pipeline.log`

## Step 4: Run Full Pipeline

Process all RNAs in the dataset:

```bash
python run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --output_dir data/rna_msa_output \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz \
    --n_cpu 8
```

**Note**: This may take several hours depending on:
- Number of RNAs in dataset
- RNA lengths
- Available computational resources
- Whether Rfam database is available

## Step 5: Verify Output

Check the generated features:

```python
import numpy as np
import json

# Load features
features = np.load('data/embeddings/rna_2d_features_rpi7317.npz')
print(f"Generated features for {len(features.files)} RNAs")

# Check a specific RNA
rna_id = features.files[0]
print(f"\nRNA: {rna_id}")
print(f"  Conservation shape: {features[rna_id]['conservation'].shape}")
print(f"  Coevolution shape: {features[rna_id]['coevolution'].shape}")
print(f"  Accessibility shape: {features[rna_id]['accessibility'].shape}")
print(f"  MSA depth: {features[rna_id]['msa_depth']}")

# Load pipeline results
with open('data/rna_msa_output/pipeline_results.json') as f:
    results = json.load(f)

# Count strategies
strategies = {}
for result in results.values():
    if result.get('success'):
        strategy = result['strategy']
        strategies[strategy] = strategies.get(strategy, 0) + 1

print("\nStrategies used:")
for strategy, count in strategies.items():
    print(f"  {strategy}: {count}")
```

## Step 6: Integrate with Model

See `examples/integrate_2d_features.py` for example of how to use the 2D features in your model.

## Troubleshooting

### Issue: "No Rfam hits found"

**Solution**: This is expected for novel RNAs. Pipeline will automatically fall back to RNAcmap-style or single-sequence prediction.

### Issue: "Low MSA depth"

**Solutions**:
1. Provide a sequence database with `--sequence_db` parameter
2. This is expected for lineage-specific RNAs
3. Single-sequence features will still be generated

### Issue: "Tool not found"

**Solution**: Install missing tools:
```bash
# Check which tools are missing
python run_rna_msa_pipeline.py --check_deps

# Install missing tools
# Ubuntu/Debian:
sudo apt-get install infernal vienna-rna mafft

# macOS:
brew install infernal viennarna mafft
```

### Issue: Pipeline is slow

**Solutions**:
1. Increase `--n_cpu` parameter
2. Run in test mode first to verify setup
3. Consider processing in batches
4. Skip Rfam search if database is not available

### Issue: Out of memory

**Solutions**:
1. Reduce `--n_cpu` parameter
2. Process RNAs in batches
3. Increase system swap space

## Expected Runtime

Approximate times (on 8-core CPU):

- **Test mode (5 RNAs)**: 2-5 minutes
- **Small dataset (100 RNAs)**: 30-60 minutes
- **RPI7317 dataset (~7000 RNAs)**: 6-12 hours

Times vary based on:
- RNA lengths
- MSA depth achieved
- Available tools
- System resources

## Next Steps

1. ✅ Verify pipeline works in test mode
2. ✅ Run on full dataset
3. ✅ Check output features
4. 📝 Integrate features into model (see `examples/integrate_2d_features.py`)
5. 📝 Train model with combined sequence + structure features
6. 📝 Compare performance with sequence-only baseline

## Getting Help

- Check `RNA_MSA_README.md` for detailed documentation
- Review log file: `rna_msa_pipeline.log`
- Check pipeline results: `data/rna_msa_output/pipeline_results.json`
