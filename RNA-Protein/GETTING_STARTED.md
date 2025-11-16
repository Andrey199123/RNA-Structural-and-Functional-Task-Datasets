# Getting Started with RNA MSA Pipeline

## Overview

This pipeline adds structural features to your RNA-protein interaction prediction model, addressing your professor's feedback about incorporating 2D features alongside sequence embeddings.

## What This Pipeline Does

1. **Generates RNA MSAs** using 4 different strategies (automatically selected)
2. **Extracts 2D features**: conservation, coevolution, and accessibility maps
3. **Outputs features** in NPZ format ready for model integration
4. **Handles all RNA types**: short, long, known families, novel sequences

## Quick Start (5 Minutes)

### Step 1: Verify Installation

```bash
cd RNA-Protein
python3 test_rna_msa_pipeline.py
```

Expected output: `✓ All tests passed! Pipeline is ready to use.`

### Step 2: Test on Small Dataset

```bash
python3 run_rna_msa_pipeline.py --test_mode --verbose
```

This processes the first 5 RNAs to verify everything works.

### Step 3: Check Output

```bash
# View results
cat data/rna_msa_output/pipeline_results.json

# Check features (if you have Python with numpy)
python3 -c "import numpy as np; f=np.load('data/embeddings/rna_2d_features_rpi7317.npz', allow_pickle=True); print(f'Generated features for {len(f.files)} RNAs')"
```

## Full Pipeline Run

### Process All RPI7317 RNAs

```bash
python3 run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --output_dir data/rna_msa_output \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz \
    --n_cpu 8
```

**Expected runtime**: 6-12 hours for ~7000 RNAs (depends on your system)

### Monitor Progress

```bash
# Watch the log file
tail -f rna_msa_pipeline.log

# Check how many RNAs processed
ls data/rna_msa_output/*_structure.npz | wc -l
```

## Understanding the Output

### 1. Feature File (NPZ)

```python
import numpy as np

# Load features
features = np.load('data/embeddings/rna_2d_features_rpi7317.npz', allow_pickle=True)

# Check an RNA
rna_id = features.files[0]
rna_features = features[rna_id].item()

print(f"RNA: {rna_id}")
print(f"  Conservation: {rna_features['conservation'].shape}")  # (length,)
print(f"  Coevolution: {rna_features['coevolution'].shape}")    # (length, length)
print(f"  Accessibility: {rna_features['accessibility'].shape}") # (length,)
print(f"  MSA depth: {rna_features['msa_depth']}")              # int
```

### 2. Pipeline Results (JSON)

```bash
# View summary
python3 -c "
import json
with open('data/rna_msa_output/pipeline_results.json') as f:
    results = json.load(f)
    
success = sum(1 for r in results.values() if r.get('success'))
print(f'Successful: {success}/{len(results)}')

strategies = {}
for r in results.values():
    if r.get('success'):
        s = r['strategy']
        strategies[s] = strategies.get(s, 0) + 1

print('\\nStrategies used:')
for s, c in strategies.items():
    print(f'  {s}: {c}')
"
```

## Integrating with Your Model

### Current Model (Sequence Only)

```python
# Your current approach
rna_embedding = rna_fm_model(rna_sequence)  # (640,)
protein_embedding = esm2_model(protein_sequence)  # (1280,)

# Predict interaction
score = model(rna_embedding, protein_embedding)
```

### Enhanced Model (Sequence + Structure)

```python
from examples.integrate_2d_features import CombinedRNAEncoder, load_features

# Initialize combined encoder
rna_encoder = CombinedRNAEncoder(
    seq_dim=640,      # RNA-FM dimension
    struct_dim=128,   # Structural feature dimension
    out_dim=128       # Output dimension
)

# Load features
features = load_features(
    rna_id,
    'data/embeddings/rna_fm_rpi7317.npz',
    'data/embeddings/rna_2d_features_rpi7317.npz'
)

# Encode RNA with both sequence and structure
rna_embedding = rna_encoder(
    features['seq_emb'],
    features['conservation'],
    features['coevolution'],
    features['accessibility']
)

# Predict interaction (same as before)
score = model(rna_embedding, protein_embedding)
```

See `examples/integrate_2d_features.py` for complete implementation.

## Visualizing Features

```bash
# Visualize features for specific RNAs
python3 examples/visualize_2d_features.py \
    --features_file data/embeddings/rna_2d_features_rpi7317.npz \
    --output_dir data/rna_msa_output/visualizations \
    --n_examples 10

# View statistics only
python3 examples/visualize_2d_features.py --statistics_only
```

This creates plots showing:
- Conservation profiles
- Accessibility profiles
- Coevolution/contact maps
- Potential binding sites
- Dataset statistics

## Common Issues & Solutions

### Issue: "No external tools found"

**Impact**: Pipeline will use fallback methods (single-sequence prediction)

**Solution** (optional, for better results):
```bash
# Ubuntu/Debian
sudo apt-get install infernal vienna-rna mafft

# macOS
brew install infernal viennarna mafft
```

**Note**: Pipeline works without these tools, just with lower MSA quality.

### Issue: "Low MSA depth"

**Impact**: Features based on single-sequence prediction

**Why**: Many lncRNAs are lineage-specific with few homologs

**Solution**: This is expected and handled automatically. Single-sequence features are still useful.

### Issue: Pipeline is slow

**Solutions**:
1. Increase `--n_cpu` (e.g., `--n_cpu 16`)
2. Run overnight or on HPC
3. Process in batches (split FASTA file)

### Issue: Out of memory

**Solutions**:
1. Reduce `--n_cpu` (e.g., `--n_cpu 2`)
2. Close other applications
3. Process in smaller batches

## Expected Results

Based on the papers (trRosetta, RNAcmap), you should see:

1. **Improved accuracy** on test sets (especially unknown RNA/protein)
2. **Better generalization** to novel sequences
3. **Domain-specific predictions** for long RNAs
4. **Interpretable features** (can visualize binding sites)

## Next Steps

### 1. Generate Features (Done ✓)
```bash
python3 run_rna_msa_pipeline.py --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta
```

### 2. Modify Your Model
- Add `CombinedRNAEncoder` from `examples/integrate_2d_features.py`
- Replace RNA encoder in your existing model
- Keep protein encoder unchanged

### 3. Train Model
```bash
# Your existing training script, but with new RNA encoder
python3 src/training/train_linkpred.py \
    --config configs/default.yaml \
    --use_2d_features  # Add this flag
```

### 4. Compare Performance
- Train baseline (sequence only)
- Train enhanced (sequence + structure)
- Compare metrics on test sets

### 5. Analyze Results
- Which features help most?
- Does it improve unknown RNA prediction?
- Visualize learned attention weights

## File Organization

```
RNA-Protein/
├── data/
│   ├── embeddings/
│   │   ├── rna_fm_rpi7317.npz              # Existing sequence embeddings
│   │   └── rna_2d_features_rpi7317.npz     # NEW: 2D structural features
│   └── rna_msa_output/                      # NEW: MSA files and results
│       ├── pipeline_results.json
│       └── visualizations/
│
├── src/rna_msa/                             # NEW: Pipeline code
├── examples/                                # NEW: Integration examples
│
├── run_rna_msa_pipeline.py                  # NEW: Main runner
├── test_rna_msa_pipeline.py                 # NEW: Test suite
│
├── RNA_MSA_README.md                        # Detailed documentation
├── QUICKSTART_RNA_MSA.md                    # Quick start guide
├── IMPLEMENTATION_SUMMARY.md                # Implementation details
└── GETTING_STARTED.md                       # This file
```

## Documentation

- **GETTING_STARTED.md** (this file): Quick start guide
- **QUICKSTART_RNA_MSA.md**: Step-by-step tutorial
- **RNA_MSA_README.md**: Comprehensive documentation
- **IMPLEMENTATION_SUMMARY.md**: Technical details

## Getting Help

1. **Check logs**: `cat rna_msa_pipeline.log`
2. **Run tests**: `python3 test_rna_msa_pipeline.py`
3. **Check results**: `cat data/rna_msa_output/pipeline_results.json`
4. **Review docs**: See documentation files above

## Timeline Estimate

- **Setup & Testing**: 30 minutes
- **Full pipeline run**: 6-12 hours
- **Model integration**: 2-4 hours
- **Training & evaluation**: 1-2 days
- **Analysis & visualization**: 1 day

**Total**: ~3-4 days to complete integration and see results

## Success Criteria

You'll know it's working when:

1. ✓ Tests pass: `python3 test_rna_msa_pipeline.py`
2. ✓ Features generated: `data/embeddings/rna_2d_features_rpi7317.npz` exists
3. ✓ Model trains: No errors with new encoder
4. ✓ Performance improves: Better metrics than baseline

## Questions?

- Review the documentation files
- Check the example scripts in `examples/`
- Look at the log file for detailed information
- Examine `pipeline_results.json` for per-RNA results

Good luck with your research! 🧬
