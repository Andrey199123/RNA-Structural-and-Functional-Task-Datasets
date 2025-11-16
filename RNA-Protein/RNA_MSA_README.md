# RNA MSA Pipeline for RNA-Protein Interaction Prediction

This pipeline generates high-quality RNA multiple sequence alignments (MSAs) and extracts 2D structural features for improving RNA-protein interaction prediction models.

## Overview

The pipeline implements four strategies for generating RNA MSAs, automatically selecting the best approach based on RNA characteristics:

1. **Infernal + Rfam**: For known/annotated RNA families
2. **RNAcmap-Style**: For novel RNAs without Rfam matches
3. **Domain Split**: For long RNAs (>1000 nt)
4. **Single-Sequence**: For orphan/low-homology RNAs

## Features Extracted

For each RNA, the pipeline generates three types of 2D features:

1. **Conservation Map** (per-position): Identifies highly conserved nucleotides important for function/structure
2. **Coevolution/Contact Map** (pairwise): Detects covarying nucleotide pairs indicating base pairing or 3D proximity
3. **Pairing/Accessibility Map**: Distinguishes paired (buried) vs unpaired (exposed) regions

These features can be used with vision-based architectures (e.g., Vision Mamba) to improve model performance.

## Installation

### Required Dependencies

```bash
# Python packages
pip install numpy scipy biopython

# Optional: ViennaRNA Python bindings for better structure prediction
pip install ViennaRNA
```

### External Tools (Optional but Recommended)

1. **Infernal** (for Rfam search):
```bash
# Ubuntu/Debian
sudo apt-get install infernal

# macOS
brew install infernal

# Or download from: http://eddylab.org/infernal/
```

2. **ViennaRNA** (for structure prediction):
```bash
# Ubuntu/Debian
sudo apt-get install vienna-rna

# macOS
brew install viennarna

# Or download from: https://www.tbi.univie.ac.at/RNA/
```

3. **MAFFT** (for sequence alignment):
```bash
# Ubuntu/Debian
sudo apt-get install mafft

# macOS
brew install mafft
```

### Optional: Rfam Database

Download Rfam covariance models for known RNA family search:

```bash
# Download and prepare Rfam database
mkdir -p data/databases
cd data/databases
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm  # Requires Infernal
cd ../..
```

## Usage

### Basic Usage

Process RNA sequences and generate 2D features:

```bash
python run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --output_dir data/rna_msa_output \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz
```

### With Rfam Database

```bash
python run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --rfam_db data/databases/Rfam.cm \
    --output_dir data/rna_msa_output \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz \
    --n_cpu 8
```

### Test Mode

Process only first 5 sequences for testing:

```bash
python run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --test_mode \
    --verbose
```

### Check Dependencies

```bash
python run_rna_msa_pipeline.py --check_deps
```

## Pipeline Strategies

### 1. Infernal + Rfam (Known RNAs)

**When used**: RNA matches a known Rfam family

**Process**:
- Search Rfam database using Infernal cmscan
- Fetch family alignment from Rfam
- Align query RNA to family MSA
- Convert to .a3m format

**Advantages**: 
- High-quality curated alignments
- Structure-aware from covariance models
- Fast for known families

**Limitations**: 
- Only works for annotated RNAs
- Limited coverage for novel sequences

### 2. RNAcmap-Style (Novel RNAs)

**When used**: RNA has no Rfam match

**Process**:
1. Predict secondary structure (SPOT-RNA or RNAfold)
2. Build custom covariance model from predicted structure
3. Search sequence databases for homologs
4. Align homologous sequences
5. Refine structure based on covariation
6. Iterate to improve MSA depth

**Advantages**:
- Works for completely novel RNAs
- Builds structure-aware MSA from scratch
- Iterative refinement improves quality

**Limitations**:
- Requires large sequence database
- Computationally intensive
- Quality depends on initial structure prediction

### 3. Domain Split (Long RNAs)

**When used**: RNA length > 1000 nt

**Process**:
- Split RNA into overlapping domains (300 nt windows, 50 nt overlap)
- Process each domain independently using RNAcmap-style or Rfam
- Generate per-domain MSAs and features
- Combine domain features into full-length representation

**Advantages**:
- Handles modular structure of long RNAs
- Each domain gets appropriate treatment
- Identifies domain-specific protein binding

**Limitations**:
- May split functional domains
- Requires careful feature combination

### 4. Single-Sequence (Orphan RNAs)

**When used**: No MSA can be generated (no homologs found)

**Process**:
- Predict structure using SPOT-RNA or RNAfold
- Compute pairing probabilities
- Estimate accessibility from structure
- Use GC content as conservation proxy

**Advantages**:
- Always produces features (fallback)
- Fast computation
- No database required

**Limitations**:
- No evolutionary information
- Lower quality than MSA-based features
- Marked as "low-depth" for downstream models

## Output Files

### Directory Structure

```
data/rna_msa_output/
├── {rna_id}_query.fa              # Query sequence
├── {rna_id}_rfam.tblout           # Rfam search results
├── {rna_id}_aligned.sto           # Aligned MSA (Stockholm)
├── {rna_id}_final.a3m             # Final MSA (A3M format)
├── {rna_id}_structure.npz         # Structure predictions
└── pipeline_results.json          # Summary of all results
```

### Feature File Format

The output NPZ file contains 2D features for each RNA:

```python
import numpy as np

# Load features
features = np.load('data/embeddings/rna_2d_features_rpi7317.npz')

# Access features for a specific RNA
rna_id = 'ENST00000123456'
conservation = features[rna_id]['conservation']      # Shape: (length,)
coevolution = features[rna_id]['coevolution']        # Shape: (length, length)
accessibility = features[rna_id]['accessibility']    # Shape: (length,)
msa_depth = features[rna_id]['msa_depth']           # Scalar
```

## Integration with Existing Model

### Adding 2D Features to Training

```python
# Load sequence embeddings (existing)
rna_seq_emb = np.load('data/embeddings/rna_fm_rpi7317.npz')

# Load 2D structural features (new)
rna_2d_features = np.load('data/embeddings/rna_2d_features_rpi7317.npz')

# Combine in model
class RNAEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.seq_encoder = ...  # Existing sequence encoder
        self.struct_encoder = ...  # New: 2D feature encoder (CNN/Vision Mamba)
        self.fusion = ...  # Combine sequence + structure
    
    def forward(self, rna_id):
        # Sequence features
        seq_feat = self.seq_encoder(rna_seq_emb[rna_id])
        
        # Structural features
        conservation = rna_2d_features[rna_id]['conservation']
        coevolution = rna_2d_features[rna_id]['coevolution']
        accessibility = rna_2d_features[rna_id]['accessibility']
        
        # Stack as image-like input (3 channels)
        struct_input = np.stack([
            conservation,
            coevolution.mean(axis=1),  # Aggregate contact map
            accessibility
        ])
        
        struct_feat = self.struct_encoder(struct_input)
        
        # Fuse features
        combined = self.fusion(seq_feat, struct_feat)
        return combined
```

## Troubleshooting

### Low MSA Depth

If many RNAs have low MSA depth (< 10 sequences):

1. **Provide sequence database**: Use `--sequence_db` with large RNA database (e.g., RNAcentral, Rfam full)
2. **Adjust E-value threshold**: Modify `evalue_threshold` in code to be more permissive
3. **Use broader search**: Include related species/transcriptomes

### Missing Dependencies

If external tools are not available:

- Pipeline will automatically fall back to single-sequence prediction
- Install tools for better results (see Installation section)
- Check installation: `python run_rna_msa_pipeline.py --check_deps`

### Memory Issues

For large datasets:

- Process in batches
- Reduce `n_cpu` parameter
- Use `--test_mode` to verify pipeline works first

### Long Runtime

To speed up:

- Increase `--n_cpu` for parallel processing
- Reduce number of refinement iterations in code
- Use pre-computed Rfam database
- Skip homology search for very long RNAs

## References

### Papers

1. **trRosetta**: Yang et al. (2020) "Improved protein structure prediction using predicted interresidue orientations"
2. **trRosettaRNA**: Baek et al. (2024) "Accurate prediction of RNA structure using evolutionary information"
3. **RNAcmap**: Zhang et al. (2021) "Accurate RNA contact prediction using evolutionary couplings"
4. **Infernal**: Nawrocki & Eddy (2013) "Infernal 1.1: 100-fold faster RNA homology searches"
5. **Rfam**: Kalvari et al. (2021) "Rfam 14: expanded coverage of metagenomic, viral and microRNA families"

### Tools

- **Infernal**: http://eddylab.org/infernal/
- **ViennaRNA**: https://www.tbi.univie.ac.at/RNA/
- **MAFFT**: https://mafft.cbrc.jp/alignment/software/
- **Rfam**: https://rfam.org/
- **SPOT-RNA**: https://sparks-lab.org/server/spot-rna/

## Next Steps

1. **Run pipeline on RPI7317 dataset**
2. **Integrate 2D features into model architecture**
3. **Experiment with vision-based encoders** (CNN, Vision Transformer, Vision Mamba)
4. **Compare performance** with sequence-only baseline
5. **Analyze which features contribute most** to improved predictions

## Contact

For questions or issues with the pipeline, please refer to the project documentation or contact the research team.
