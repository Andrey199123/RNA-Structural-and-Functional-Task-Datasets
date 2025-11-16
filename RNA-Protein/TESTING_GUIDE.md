# Testing Guide for RNA MSA Pipeline

## Two Types of Tests

### 1. Functionality Tests (`test_rna_msa_pipeline.py`)

**What it tests**: Code works correctly
- ✓ Dependencies are available
- ✓ RNA sequences load properly
- ✓ Structure prediction runs
- ✓ Features are extracted
- ✓ Pipeline completes without errors
- ✓ Files are saved correctly

**Purpose**: Verify installation and basic functionality

**Run**: `python3 test_rna_msa_pipeline.py`

**Result**: All tests pass ✅ (already verified)

### 2. Quality Tests (`test_msa_quality.py` & `compare_strategies.py`)

**What it tests**: MSA quality is good
- MSA depth (number of sequences)
- Sequence diversity (not all identical)
- Conservation patterns
- Coevolution signals
- Strategy comparison

**Purpose**: Determine which approach gives best results

**Requires**: Databases and external tools for meaningful results

## Why You're Right

You're absolutely correct that the current test suite only checks if the code *runs*, not if it produces *high-quality* MSAs. Here's why:

### Current Situation

**Without databases**:
- Pipeline uses single-sequence fallback
- MSA depth = 0 (no actual MSA)
- Quality = "fair" (structure prediction only)
- Tests pass ✅ but quality is limited

**With databases**:
- Pipeline can use Infernal+Rfam and RNAcmap strategies
- MSA depth = 10-1000+ sequences
- Quality = "good" to "excellent"
- Much better evolutionary information

## Recommended Testing Workflow

### Phase 1: Verify Installation (Done ✅)

```bash
# Check code works
python3 test_rna_msa_pipeline.py

# Expected: All 6 tests pass
```

**Status**: ✅ Complete - code is functional

### Phase 2: Setup Databases (Do This Next)

```bash
# Install external tools
sudo apt-get install infernal vienna-rna mafft  # Ubuntu
# or
brew install infernal viennarna mafft  # macOS

# Download Rfam database (~2 GB)
mkdir -p data/databases
cd data/databases
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
cd ../..
```

**Time**: 10-20 minutes

**See**: `SETUP_DATABASES.md` for detailed instructions

### Phase 3: Compare Strategies

```bash
# Test different strategies on sample RNAs
python3 compare_strategies.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --n_test 20 \
    --rfam_db data/databases/Rfam.cm \
    --output_dir data/strategy_comparison
```

**Output**: Shows which strategy works best for your RNAs

**Example output**:
```
Strategy Comparison Summary
====================================

single_sequence:
  Success rate: 20/20 (100.0%)
  Average MSA depth: 0.0
  Average quality score: 0.20

infernal_rfam:
  Success rate: 5/20 (25.0%)
  Average MSA depth: 145.2
  Average quality score: 0.85

Recommendation: Use automatic strategy selection
(pipeline chooses best strategy per RNA)
```

### Phase 4: Quality Assessment

```bash
# Run quality tests
python3 test_msa_quality.py
```

**Output**: Detailed quality metrics for generated MSAs

### Phase 5: Full Pipeline with Quality Monitoring

```bash
# Run full pipeline with databases
python3 run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --rfam_db data/databases/Rfam.cm \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz \
    --n_cpu 8

# Check results
python3 -c "
import json
with open('data/rna_msa_output/pipeline_results.json') as f:
    results = json.load(f)

# Count by strategy
strategies = {}
depths = []
for r in results.values():
    if r.get('success'):
        s = r['strategy']
        strategies[s] = strategies.get(s, 0) + 1
        depths.append(r.get('msa_depth', 0))

print('Strategies used:')
for s, c in strategies.items():
    print(f'  {s}: {c}')

import numpy as np
print(f'\nMSA depth statistics:')
print(f'  Mean: {np.mean(depths):.1f}')
print(f'  Median: {np.median(depths):.1f}')
print(f'  Max: {max(depths)}')
print(f'  With MSA (depth>0): {sum(1 for d in depths if d > 0)}')
"
```

## Understanding Quality Metrics

### MSA Depth
- **0**: Single sequence (no MSA)
- **1-5**: Very shallow (limited information)
- **5-20**: Shallow (some evolutionary signal)
- **20-50**: Moderate (good signal)
- **50+**: Deep (excellent signal)
- **100+**: Very deep (optimal)

### Quality Ratings
- **Excellent** (score ≥ 0.8): Deep MSA, good diversity, clear signals
- **Good** (score ≥ 0.6): Moderate MSA, reasonable quality
- **Fair** (score ≥ 0.4): Shallow MSA or single-sequence with structure
- **Poor** (score ≥ 0.2): Very shallow MSA
- **Very Poor** (score < 0.2): No useful information

### Expected Results for lncRNAs

**With Rfam only**:
- ~20-30% match known families → excellent quality
- ~70-80% use single-sequence → fair quality
- Overall: Mixed quality

**With Rfam + sequence database**:
- ~20-30% match known families → excellent quality
- ~30-50% find homologs → good to excellent quality
- ~20-40% use single-sequence → fair quality
- Overall: Good quality

## What to Test

### Minimum (No Databases)
```bash
# Just verify code works
python3 test_rna_msa_pipeline.py
```
**Result**: Code functional, but quality limited

### Recommended (With Rfam)
```bash
# 1. Setup
# Download Rfam database (see SETUP_DATABASES.md)

# 2. Compare strategies
python3 compare_strategies.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --n_test 20 \
    --rfam_db data/databases/Rfam.cm

# 3. Run full pipeline
python3 run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --rfam_db data/databases/Rfam.cm \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz
```
**Result**: Good quality for known families, fair for novel

### Optimal (With Both Databases)
```bash
# Add sequence database for homology search
python3 run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --rfam_db data/databases/Rfam.cm \
    --sequence_db data/databases/rfam_sequences.fasta \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz
```
**Result**: Good to excellent quality for most RNAs

## Your Question: Which Strategy is Best?

**Answer**: It depends on the RNA!

- **Known families** (tRNA, rRNA, etc.): Infernal+Rfam is best
- **Novel with homologs**: RNAcmap-style is best
- **Long RNAs**: Domain split is best
- **Orphan RNAs**: Single-sequence is only option

**Solution**: Let the pipeline choose automatically!

The pipeline's `determine_strategy()` function automatically selects the best approach for each RNA based on:
- RNA length
- Rfam matches
- Homolog availability

You don't need to choose manually.

## Practical Recommendation

### For Your Research

1. **Start simple** (today):
   ```bash
   # Verify code works
   python3 test_rna_msa_pipeline.py
   ```

2. **Add Rfam** (this week):
   ```bash
   # Download Rfam (~20 min)
   # See SETUP_DATABASES.md
   
   # Test with Rfam
   python3 compare_strategies.py --rfam_db data/databases/Rfam.cm
   ```

3. **Run full pipeline** (overnight):
   ```bash
   python3 run_rna_msa_pipeline.py \
       --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
       --rfam_db data/databases/Rfam.cm
   ```

4. **Evaluate quality**:
   ```bash
   # Check MSA depths in results
   # Visualize features
   python3 examples/visualize_2d_features.py
   ```

5. **Integrate and train**:
   ```bash
   # Use features in your model
   # Compare with baseline
   ```

### Expected Timeline

- **Setup databases**: 30 minutes
- **Strategy comparison**: 10 minutes
- **Full pipeline**: 6-12 hours (overnight)
- **Quality evaluation**: 30 minutes
- **Model integration**: 2-4 hours

**Total**: ~1-2 days to high-quality features

## Summary

**Current status**: ✅ Code works (functionality tests pass)

**Next step**: Download Rfam database for quality improvement

**Goal**: Generate high-quality MSAs with evolutionary information

**How**: Follow `SETUP_DATABASES.md` → Run `compare_strategies.py` → Run full pipeline

You're absolutely right to focus on quality, not just functionality!
