# RNA MSA Pipeline - Implementation Status

## Overview

This document describes the current implementation status of the RNA MSA pipeline for extracting 2D structural features from RNA sequences.

## What's Implemented ✅

### 1. Real Rfam Integration
- ✅ Extract family alignments from Rfam.seed file
- ✅ Align query RNA to family using cmalign
- ✅ Fallback to mock alignments if Rfam.seed not available
- ✅ Parse Stockholm format from Rfam
- ✅ Extract real evolutionary signal from family MSAs

**Files:**
- `src/rna_msa/infernal_rfam.py` - Functions: `extract_family_from_seed()`, `align_to_family_with_cmalign()`

### 2. MSA Cleaning & Quality Control
- ✅ Remove near-duplicate sequences (>95% identity)
- ✅ Filter partial sequences (low coverage)
- ✅ Remove sequences with too many gaps
- ✅ Compute MSA statistics (depth, identity, coverage)
- ✅ Re-alignment with MAFFT (structure-aware)

**Files:**
- `src/rna_msa/msa_utils.py` - Functions: `clean_msa()`, `compute_msa_statistics()`

### 3. Format Conversion
- ✅ Convert Stockholm → A3M format
- ✅ Proper A3M formatting:
  - Query sequence first (no gaps)
  - Other sequences: '-' for deletions, lowercase for insertions
- ✅ Validate format for downstream models

**Files:**
- `src/rna_msa/msa_utils.py` - Function: `convert_stockholm_to_a3m()`

### 4. 2D Feature Extraction
- ✅ Conservation map (per-position)
- ✅ Coevolution/contact map (pairwise, using mutual information + APC)
- ✅ Accessibility/pairing map
- ✅ Support for multiple MSA formats (Stockholm, A3M, FASTA)

**Files:**
- `src/rna_msa/feature_extraction.py`

### 5. Testing Infrastructure
- ✅ Test script for RPI7317 dataset
- ✅ Dataset characterization (length distribution, RNA types)
- ✅ Success rate measurement per strategy
- ✅ MSA depth and quality benchmarking

**Files:**
- `test_rpi7317_dataset.py`

## What's Partially Implemented ⚠️

### 1. RNAcmap-Style Pipeline
- ✅ Basic structure prediction with SPOT-RNA
- ✅ Custom CM building from predicted structure
- ⚠️ Iterative refinement (needs sequence database)
- ⚠️ Structure update based on covariation (needs testing)

**Status:** Works for single-sequence prediction, but iterative refinement requires sequence database.

### 2. Domain Splitting for Long RNAs
- ✅ Simple sliding window approach
- ⚠️ Structural domain prediction (needs RNAStructure/ViennaRNA integration)
- ⚠️ Intelligent domain boundary detection
- ✅ Feature combination from multiple domains

**Status:** Basic implementation works, but could be improved with better domain detection.

### 3. Sequence Database
- ⚠️ Not included in repository (too large)
- ✅ Documentation for downloading RNAcentral
- ⚠️ Formatting and indexing scripts needed

**Status:** Users need to download and set up their own sequence database.

## Usage

### Quick Test on tRNA
```bash
conda activate rna_msa

python3 run_rna_msa_pipeline.py \
  --rna_fasta data/test_known_rnas.fasta \
  --rfam_db data/databases/Rfam.cm \
  --rfam_seed data/databases/Rfam.seed \
  --output_features data/embeddings/test_features.npz \
  --verbose
```

### Test on RPI7317 Dataset
```bash
# Analyze dataset only
python3 test_rpi7317_dataset.py \
  --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
  --analyze_only

# Test on first 10 RNAs
python3 test_rpi7317_dataset.py \
  --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
  --rfam_db data/databases/Rfam.cm \
  --rfam_seed data/databases/Rfam.seed \
  --max_test 10

# Full dataset test
python3 test_rpi7317_dataset.py \
  --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
  --rfam_db data/databases/Rfam.cm \
  --rfam_seed data/databases/Rfam.seed \
  --output_dir data/rna_msa_test
```

## Requirements

### Required Tools
- **Infernal** (cmscan, cmalign, cmbuild) - For Rfam search and alignment
- **MAFFT** - For MSA re-alignment
- **ViennaRNA** (RNAfold) - For structure prediction
- **Python packages**: numpy, biopython, scipy

### Optional Tools
- **SPOT-RNA** - For better structure prediction
- **RNAStructure** - For domain detection
- **R-Coffee/LocARNA** - For structure-aware alignment

### Database Files
- **Rfam.cm** (~311 MB) - Covariance models
- **Rfam.seed** (~1.4 GB) - Family seed alignments
- **Sequence database** (optional, for RNAcmap) - RNAcentral or custom

Download instructions in `data/databases/README.md`

## Current Limitations

1. **Rfam.seed not included** - Too large for GitHub, users must download
2. **Mock alignments as fallback** - If Rfam.seed not available, uses mock data
3. **No sequence database** - RNAcmap iterative refinement needs external database
4. **Simple domain splitting** - Uses fixed windows instead of structural domains
5. **Limited testing** - Tested on tRNA, needs full RPI7317 validation

## Next Steps

### To Satisfy Professor's Email (Priority)
1. ✅ Implement real Rfam integration - **DONE**
2. ✅ Add MSA cleaning and quality control - **DONE**
3. ✅ Implement .a3m format conversion - **DONE**
4. ⏳ Test on RPI7317 dataset - **READY TO RUN**
5. ⏳ Document results and challenges - **AFTER TESTING**

### Future Improvements
1. Download and integrate sequence database for RNAcmap
2. Improve domain detection for long RNAs
3. Add more coevolution methods (DCA, plmc)
4. Integrate with protein 2D features
5. Train vision-based model (Vision Mamba)

## Testing Checklist

- [x] Test on short RNA (tRNA, 76nt)
- [ ] Test on RPI7317 dataset
- [ ] Characterize RNA length distribution
- [ ] Measure success rate per strategy
- [ ] Benchmark MSA depth and quality
- [ ] Identify problematic cases (long RNAs, orphans)
- [ ] Document gap problem with long RNAs
- [ ] Generate report for professor

## Files Modified/Created

### Core Implementation
- `src/rna_msa/infernal_rfam.py` - Real Rfam integration
- `src/rna_msa/msa_utils.py` - MSA cleaning and .a3m conversion
- `src/rna_msa/pipeline.py` - Added rfam_seed_path parameter
- `src/rna_msa/feature_extraction.py` - Path handling fix

### Testing & Documentation
- `test_rpi7317_dataset.py` - Dataset testing script
- `run_rna_msa_pipeline.py` - Updated with rfam_seed parameter
- `RNA_MSA_IMPLEMENTATION_STATUS.md` - This file
- `data/databases/README.md` - Database download instructions

## Contact

For questions about this implementation, contact the research team.

Last updated: November 2025
