# RNA MSA Pipeline - Documentation Index

## 🚀 Start Here

**New to this pipeline?** Start with [GETTING_STARTED.md](GETTING_STARTED.md) for a 5-minute quick start.

## 📚 Documentation Guide

### For Quick Start
- **[GETTING_STARTED.md](GETTING_STARTED.md)** - 5-minute quick start guide
  - Verify installation
  - Run test mode
  - Check output
  - Integrate with model

### For Step-by-Step Tutorial
- **[QUICKSTART_RNA_MSA.md](QUICKSTART_RNA_MSA.md)** - Detailed tutorial
  - Installation instructions
  - Troubleshooting guide
  - Expected runtimes
  - Next steps

### For Complete Documentation
- **[RNA_MSA_README.md](RNA_MSA_README.md)** - Comprehensive documentation
  - Full pipeline description
  - All four strategies explained
  - Feature extraction details
  - Integration examples
  - References and papers

### For Technical Details
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** - Implementation details
  - Project structure
  - File descriptions
  - Function documentation
  - Output formats
  - API reference

### For Project Overview
- **[PROJECT_SUMMARY.md](PROJECT_SUMMARY.md)** - High-level summary
  - What was implemented
  - Why it matters
  - How to use it
  - Expected benefits
  - Timeline to results

## 🎯 Quick Links by Task

### I want to...

**...understand what this does**
→ Read [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md)

**...get started quickly**
→ Follow [GETTING_STARTED.md](GETTING_STARTED.md)

**...learn step-by-step**
→ Follow [QUICKSTART_RNA_MSA.md](QUICKSTART_RNA_MSA.md)

**...understand the science**
→ Read [RNA_MSA_README.md](RNA_MSA_README.md)

**...see the code structure**
→ Read [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)

**...integrate with my model**
→ See `examples/integrate_2d_features.py`

**...visualize the features**
→ Run `examples/visualize_2d_features.py`

**...troubleshoot issues**
→ Check [QUICKSTART_RNA_MSA.md](QUICKSTART_RNA_MSA.md) Troubleshooting section

## 📁 File Structure

```
RNA-Protein/
│
├── Documentation (6 files)
│   ├── README_RNA_MSA.md              ← You are here
│   ├── GETTING_STARTED.md             ← Start here!
│   ├── QUICKSTART_RNA_MSA.md          ← Tutorial
│   ├── RNA_MSA_README.md              ← Full docs
│   ├── IMPLEMENTATION_SUMMARY.md      ← Technical
│   └── PROJECT_SUMMARY.md             ← Overview
│
├── Pipeline Code
│   └── src/rna_msa/                   ← 7 Python modules
│       ├── pipeline.py                ← Main orchestrator
│       ├── infernal_rfam.py          ← Rfam strategy
│       ├── rnacmap_style.py          ← Novel RNA strategy
│       ├── domain_split.py           ← Long RNA strategy
│       ├── single_sequence.py        ← Fallback strategy
│       ├── feature_extraction.py     ← 2D features
│       └── msa_utils.py              ← Utilities
│
├── Examples
│   ├── integrate_2d_features.py      ← Model integration
│   └── visualize_2d_features.py      ← Visualization
│
├── Scripts
│   ├── run_rna_msa_pipeline.py       ← Main runner
│   └── test_rna_msa_pipeline.py      ← Test suite
│
└── Dependencies
    └── requirements_rna_msa.txt       ← Python packages
```

## 🔬 What This Pipeline Does

1. **Generates RNA MSAs** using 4 strategies (auto-selected)
2. **Extracts 2D features**: conservation, coevolution, accessibility
3. **Outputs features** ready for model integration
4. **Handles all RNA types**: short, long, known, novel, orphan

## 🎓 Scientific Background

This pipeline implements methods from:
- **trRosetta** (Yang et al., 2020) - Protein structure prediction
- **trRosettaRNA** (Baek et al., 2024) - RNA structure prediction
- **RNAcmap** (Zhang et al., 2021) - RNA contact prediction
- **Infernal** (Nawrocki & Eddy, 2013) - RNA homology search
- **Rfam** (Kalvari et al., 2021) - RNA family database

## ⚡ Quick Commands

```bash
# Test installation
python3 test_rna_msa_pipeline.py

# Test mode (first 5 RNAs)
python3 run_rna_msa_pipeline.py --test_mode

# Full run
python3 run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz

# Visualize
python3 examples/visualize_2d_features.py
```

## 📊 Expected Results

- **MSA depth**: Varies by RNA (0-1000+ sequences)
- **Feature quality**: High for known RNAs, moderate for novel
- **Runtime**: 6-12 hours for ~7000 RNAs
- **Success rate**: ~95%+ (fallback ensures all RNAs get features)

## 🎯 Integration with Your Model

### Current (Sequence Only)
```python
rna_emb = rna_fm_model(sequence)  # (640,)
```

### Enhanced (Sequence + Structure)
```python
from examples.integrate_2d_features import CombinedRNAEncoder

encoder = CombinedRNAEncoder(seq_dim=640, struct_dim=128)
rna_emb = encoder(seq_emb, conservation, coevolution, accessibility)
```

See `examples/integrate_2d_features.py` for complete code.

## ✅ Testing Status

All tests pass:
- ✓ Dependencies check
- ✓ RNA sequence loading
- ✓ Structure prediction
- ✓ Feature extraction
- ✓ Pipeline integration
- ✓ Feature saving/loading

Run: `python3 test_rna_msa_pipeline.py`

## 🔧 Dependencies

### Required
- Python 3.7+
- numpy, scipy, biopython

### Optional (for better results)
- Infernal (RNA homology search)
- ViennaRNA (structure prediction)
- MAFFT (sequence alignment)

**Note**: Pipeline works without optional tools using fallback methods.

## 📈 Timeline

- **Setup**: 30 minutes
- **Pipeline run**: 6-12 hours (overnight)
- **Integration**: 2-4 hours
- **Training**: 1-2 days
- **Analysis**: 1 day

**Total**: ~3-4 days to results

## 🆘 Getting Help

1. **Check documentation** - See files above
2. **Run tests** - `python3 test_rna_msa_pipeline.py`
3. **Check logs** - `cat rna_msa_pipeline.log`
4. **View results** - `cat data/rna_msa_output/pipeline_results.json`

## 🎉 Ready to Start?

1. Read [GETTING_STARTED.md](GETTING_STARTED.md)
2. Run `python3 test_rna_msa_pipeline.py`
3. Try `python3 run_rna_msa_pipeline.py --test_mode`
4. Run full pipeline
5. Integrate into your model

**Good luck with your research!** 🧬🔬

---

*This pipeline was implemented to address your professor's feedback about adding structural features to your RNA-protein interaction prediction model.*
