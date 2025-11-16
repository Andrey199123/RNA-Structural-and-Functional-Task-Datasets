# RNA MSA Pipeline - Implementation Summary

## What Has Been Implemented

I've implemented a comprehensive RNA Multiple Sequence Alignment (MSA) pipeline that addresses your professor's requirements for adding structural features to your RNA-protein interaction prediction model.

## Project Structure

```
RNA-Protein/
├── src/rna_msa/                          # Main pipeline package
│   ├── __init__.py                       # Package initialization
│   ├── pipeline.py                       # Main pipeline orchestrator
│   ├── infernal_rfam.py                  # Strategy 1: Infernal + Rfam
│   ├── rnacmap_style.py                  # Strategy 2: RNAcmap-style
│   ├── domain_split.py                   # Strategy 3: Domain splitting
│   ├── single_sequence.py                # Strategy 4: Single-sequence fallback
│   ├── feature_extraction.py             # 2D feature extraction
│   └── msa_utils.py                      # Utility functions
│
├── examples/                              # Example scripts
│   ├── integrate_2d_features.py          # Model integration example
│   └── visualize_2d_features.py          # Feature visualization
│
├── run_rna_msa_pipeline.py               # Main runner script
├── test_rna_msa_pipeline.py              # Test suite
│
├── RNA_MSA_README.md                     # Comprehensive documentation
├── QUICKSTART_RNA_MSA.md                 # Quick start guide
├── requirements_rna_msa.txt              # Python dependencies
└── IMPLEMENTATION_SUMMARY.md             # This file
```

## Four Strategies Implemented

### 1. Infernal + Rfam (Known RNAs)
- **File**: `src/rna_msa/infernal_rfam.py`
- **When used**: RNA matches a known Rfam family
- **Process**: 
  - Search Rfam database using Infernal cmscan
  - Fetch family alignment
  - Align query to family MSA
  - Convert to .a3m format
- **Key functions**:
  - `search_rfam()`: Search Rfam database
  - `fetch_rfam_alignment()`: Get family alignment
  - `convert_stockholm_to_a3m()`: Format conversion

### 2. RNAcmap-Style (Novel RNAs)
- **File**: `src/rna_msa/rnacmap_style.py`
- **When used**: RNA has no Rfam match
- **Process**:
  - Predict secondary structure (SPOT-RNA/RNAfold)
  - Build custom covariance model
  - Search for homologs
  - Iteratively refine MSA
- **Key functions**:
  - `predict_secondary_structure()`: Structure prediction
  - `build_custom_cm()`: Build covariance model
  - `iterative_msa_refinement()`: Refine MSA iteratively

### 3. Domain Split (Long RNAs)
- **File**: `src/rna_msa/domain_split.py`
- **When used**: RNA length > 1000 nt
- **Process**:
  - Split into overlapping domains (300 nt windows)
  - Process each domain independently
  - Combine domain features
- **Key functions**:
  - `split_into_domains()`: Split RNA into domains
  - `predict_domain_boundaries()`: Smart domain detection
  - `merge_overlapping_domains()`: Handle overlaps

### 4. Single-Sequence (Orphan RNAs)
- **File**: `src/rna_msa/single_sequence.py`
- **When used**: No MSA can be generated
- **Process**:
  - Predict structure (SPOT-RNA/RNAfold)
  - Compute pairing probabilities
  - Estimate accessibility
  - Use GC content as conservation proxy
- **Key functions**:
  - `predict_structure_spotRNA()`: Structure prediction
  - `compute_accessibility()`: Accessibility from structure
  - `predict_binding_sites()`: Identify potential binding sites

## Feature Extraction

### File: `src/rna_msa/feature_extraction.py`

Extracts three types of 2D features:

1. **Conservation Map** (per-position)
   - Shannon entropy-based conservation
   - Identifies functionally important positions
   - Shape: (length,)

2. **Coevolution/Contact Map** (pairwise)
   - Mutual information between positions
   - APC correction for background
   - Identifies base pairing and 3D contacts
   - Shape: (length, length)

3. **Pairing/Accessibility Map** (per-position)
   - Distinguishes paired vs unpaired regions
   - Identifies potential protein binding sites
   - Shape: (length,)

**Key functions**:
- `extract_2d_features()`: Main extraction function
- `compute_conservation_map()`: Conservation scores
- `compute_coevolution_map()`: Coevolution matrix
- `compute_accessibility_map()`: Accessibility scores

## Main Pipeline

### File: `src/rna_msa/pipeline.py`

**Class: `RNAMSAPipeline`**

Main orchestrator that:
1. Automatically selects appropriate strategy
2. Processes individual RNAs or entire datasets
3. Handles fallbacks when strategies fail
4. Extracts and saves 2D features

**Key methods**:
- `determine_strategy()`: Auto-select strategy
- `process_rna()`: Process single RNA
- `process_dataset()`: Process entire dataset
- `_extract_all_features()`: Extract features for all RNAs

## Utility Functions

### File: `src/rna_msa/msa_utils.py`

- `clean_msa()`: Remove redundant/low-quality sequences
- `realign_msa()`: Structure-aware realignment
- `compute_msa_statistics()`: MSA quality metrics
- `check_dependencies()`: Verify tool installation

## Usage Scripts

### 1. Main Runner: `run_rna_msa_pipeline.py`

Command-line interface for the pipeline:

```bash
# Basic usage
python run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --output_dir data/rna_msa_output \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz

# With Rfam database
python run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --rfam_db data/databases/Rfam.cm \
    --n_cpu 8

# Test mode (first 5 sequences)
python run_rna_msa_pipeline.py --test_mode --verbose
```

### 2. Test Suite: `test_rna_msa_pipeline.py`

Comprehensive tests to verify installation:

```bash
python test_rna_msa_pipeline.py
```

Tests:
- Dependency checking
- RNA sequence loading
- Structure prediction
- Feature extraction
- Pipeline integration
- Feature saving/loading

### 3. Model Integration: `examples/integrate_2d_features.py`

Shows how to use 2D features in your model:

```python
# Load features
features = load_features(rna_id, seq_emb_file, struct_feat_file)

# Use combined encoder
model = CombinedRNAEncoder(seq_dim=640, struct_dim=128)
output = model(
    features['seq_emb'],
    features['conservation'],
    features['coevolution'],
    features['accessibility']
)
```

### 4. Visualization: `examples/visualize_2d_features.py`

Visualize extracted features:

```bash
# Visualize specific RNAs
python examples/visualize_2d_features.py \
    --features_file data/embeddings/rna_2d_features_rpi7317.npz \
    --output_dir data/rna_msa_output/visualizations

# Statistics only
python examples/visualize_2d_features.py --statistics_only
```

## Output Format

### Feature NPZ File

```python
import numpy as np

features = np.load('data/embeddings/rna_2d_features_rpi7317.npz')

# For each RNA:
rna_features = features[rna_id].item()
# {
#     'conservation': array (length,),
#     'coevolution': array (length, length),
#     'accessibility': array (length,),
#     'msa_depth': int
# }
```

### Pipeline Results JSON

```json
{
  "rna_id": {
    "success": true,
    "strategy": "rnacmap_style",
    "length": 1234,
    "msa_depth": 45,
    "msa_path": "path/to/msa.a3m"
  }
}
```

## Dependencies

### Required (Python)
- numpy >= 1.20.0
- scipy >= 1.7.0
- biopython >= 1.79

### Optional (External Tools)
- **Infernal**: For Rfam search (cmscan, cmbuild, cmsearch)
- **ViennaRNA**: For structure prediction (RNAfold)
- **MAFFT**: For sequence alignment

### Installation

```bash
# Python packages
pip install -r requirements_rna_msa.txt

# External tools (Ubuntu/Debian)
sudo apt-get install infernal vienna-rna mafft

# External tools (macOS)
brew install infernal viennarna mafft
```

## Key Features

1. **Automatic Strategy Selection**: Pipeline automatically chooses the best approach based on RNA characteristics

2. **Robust Fallbacks**: If a strategy fails, automatically falls back to simpler methods

3. **Handles All RNA Types**:
   - Short RNAs (< 100 nt)
   - Medium RNAs (100-1000 nt)
   - Long RNAs (> 1000 nt)
   - Known families (Rfam)
   - Novel sequences
   - Orphan RNAs

4. **Quality Control**:
   - MSA cleaning (remove duplicates, low-quality sequences)
   - Gap filtering
   - Coverage checking

5. **Comprehensive Output**:
   - MSA files (.a3m, Stockholm)
   - 2D feature maps (NPZ)
   - Pipeline results (JSON)
   - Detailed logs

## Integration with Existing Model

The 2D features can be integrated into your existing model in several ways:

### Option 1: Feature Concatenation
Concatenate 2D features with sequence embeddings

### Option 2: Dual Encoder
- Sequence encoder for RNA-FM embeddings
- Structure encoder for 2D features (CNN/Vision Mamba)
- Fusion layer to combine

### Option 3: Multi-Modal Attention
Use attention mechanism to fuse sequence and structure

See `examples/integrate_2d_features.py` for implementation examples.

## Next Steps

1. **Install dependencies**:
   ```bash
   pip install -r requirements_rna_msa.txt
   ```

2. **Run tests**:
   ```bash
   python test_rna_msa_pipeline.py
   ```

3. **Test on small dataset**:
   ```bash
   python run_rna_msa_pipeline.py --test_mode --verbose
   ```

4. **Run on full RPI7317 dataset**:
   ```bash
   python run_rna_msa_pipeline.py \
       --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
       --output_features data/embeddings/rna_2d_features_rpi7317.npz \
       --n_cpu 8
   ```

5. **Visualize features**:
   ```bash
   python examples/visualize_2d_features.py
   ```

6. **Integrate into model**:
   - Modify your model architecture to use `CombinedRNAEncoder`
   - Train with both sequence and structure features
   - Compare performance with sequence-only baseline

## Expected Performance Improvements

Based on the trRosetta and RNAcmap papers, adding 2D structural features should provide:

- **Better generalization**: Structure is more conserved than sequence
- **Improved unknown RNA prediction**: Structural features help with novel sequences
- **Domain-specific binding**: Long RNA domains can be analyzed separately
- **Vision model compatibility**: 2D features work with Vision Mamba, CNNs, etc.

## Documentation

- **RNA_MSA_README.md**: Comprehensive documentation
- **QUICKSTART_RNA_MSA.md**: Quick start guide
- **IMPLEMENTATION_SUMMARY.md**: This file
- Code comments: Extensive inline documentation

## Contact & Support

For questions or issues:
1. Check the documentation files
2. Review the log file: `rna_msa_pipeline.log`
3. Check pipeline results: `data/rna_msa_output/pipeline_results.json`
4. Run tests: `python test_rna_msa_pipeline.py`

## References

Implementation based on:
- trRosetta (Yang et al., 2020)
- trRosettaRNA (Baek et al., 2024)
- RNAcmap (Zhang et al., 2021)
- Infernal (Nawrocki & Eddy, 2013)
- Rfam (Kalvari et al., 2021)
