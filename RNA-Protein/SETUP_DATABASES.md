# Setting Up Databases for RNA MSA Pipeline

## Overview

For **high-quality MSA generation**, you need:

1. **Rfam database** - For known RNA families (Strategy 1: Infernal+Rfam)
2. **Sequence database** - For homology search (Strategy 2: RNAcmap-style)
3. **External tools** - Infernal, ViennaRNA, MAFFT

**Without these**: Pipeline uses single-sequence fallback (lower quality but always works)

**With these**: Pipeline generates high-quality MSAs with evolutionary information

## Priority Setup

### Essential (Do This First)

**1. Install External Tools**

These are needed for all strategies:

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y infernal vienna-rna mafft

# macOS
brew install infernal viennarna mafft

# Verify installation
cmscan -h
RNAfold --version
mafft --version
```

**2. Download Rfam Database**

This enables Strategy 1 (Infernal+Rfam) for known RNA families:

```bash
# Create database directory
mkdir -p RNA-Protein/data/databases
cd RNA-Protein/data/databases

# Download Rfam covariance models (~500 MB compressed, ~2 GB uncompressed)
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz

# Decompress
gunzip Rfam.cm.gz

# Press database (creates index for fast searching)
cmpress Rfam.cm

# This creates:
# - Rfam.cm.i1m
# - Rfam.cm.i1p
# - Rfam.cm.i1f
# - Rfam.cm.i1i

cd ../../..
```

**Expected time**: 10-20 minutes (depending on internet speed)

**Disk space**: ~2 GB

### Optional (For Better Results)

**3. Download Sequence Database**

This enables Strategy 2 (RNAcmap-style) for novel RNAs:

**Option A: RNAcentral (Recommended)**

```bash
cd RNA-Protein/data/databases

# Download RNAcentral sequences (~50 GB compressed)
# WARNING: This is LARGE!
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz

# Decompress
gunzip rnacentral_active.fasta.gz

cd ../../..
```

**Option B: Rfam Full Sequences (Smaller)**

```bash
cd RNA-Protein/data/databases

# Download Rfam full sequences (~5 GB)
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF*.fa.gz

# Decompress all
gunzip *.fa.gz

# Concatenate into single file
cat RF*.fa > rfam_sequences.fasta

# Clean up
rm RF*.fa

cd ../../..
```

**Option C: Custom Database**

Use your own RNA sequences:
- Transcriptome assemblies
- Genome annotations
- Public RNA databases

```bash
# Just provide path to FASTA file
# Pipeline will use it for homology search
```

## Verification

### Check Installation

```bash
cd RNA-Protein

# Check tools
python3 run_rna_msa_pipeline.py --check_deps

# Expected output:
#   ✓ infernal
#   ✓ mafft
#   ✓ RNAfold
#   ✓ python_RNA (optional)
```

### Check Databases

```bash
# Check Rfam
ls -lh data/databases/Rfam.cm*

# Should see:
# Rfam.cm      (~2 GB)
# Rfam.cm.i1m  (~100 MB)
# Rfam.cm.i1p  (~50 MB)
# Rfam.cm.i1f  (~10 MB)
# Rfam.cm.i1i  (~1 MB)

# Check sequence database (if downloaded)
ls -lh data/databases/*.fasta
```

### Test with Databases

```bash
# Test with Rfam database
python3 run_rna_msa_pipeline.py \
    --test_mode \
    --rfam_db data/databases/Rfam.cm \
    --verbose

# Test with sequence database (if available)
python3 run_rna_msa_pipeline.py \
    --test_mode \
    --rfam_db data/databases/Rfam.cm \
    --sequence_db data/databases/rfam_sequences.fasta \
    --verbose
```

## Usage with Databases

### Full Pipeline with Rfam Only

```bash
python3 run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --rfam_db data/databases/Rfam.cm \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz \
    --n_cpu 8
```

### Full Pipeline with Both Databases

```bash
python3 run_rna_msa_pipeline.py \
    --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta \
    --rfam_db data/databases/Rfam.cm \
    --sequence_db data/databases/rfam_sequences.fasta \
    --output_features data/embeddings/rna_2d_features_rpi7317.npz \
    --n_cpu 8
```

## Expected Quality Improvements

### Without Databases (Baseline)
- **Strategy used**: Single-sequence only
- **MSA depth**: 0 (no MSA)
- **Quality**: Fair (structure prediction only)
- **Success rate**: 100% (always works)

### With Rfam Database
- **Strategy used**: Infernal+Rfam for known families, single-sequence for others
- **MSA depth**: 10-1000+ for known families, 0 for novel
- **Quality**: Good for known families, fair for novel
- **Success rate**: ~100%
- **Coverage**: ~20-30% of lncRNAs match Rfam families

### With Rfam + Sequence Database
- **Strategy used**: All strategies available
- **MSA depth**: 5-1000+ depending on RNA
- **Quality**: Good to excellent for most RNAs
- **Success rate**: ~95%+
- **Coverage**: ~60-80% of RNAs get meaningful MSAs

## Disk Space Requirements

```
Minimum (no databases):
- Pipeline code: ~1 MB
- Output per RNA: ~100 KB
- Total for 7000 RNAs: ~1 GB

With Rfam:
- Rfam database: ~2 GB
- Total: ~3 GB

With Rfam + RNAcentral:
- RNAcentral: ~50 GB
- Total: ~53 GB

With Rfam + Rfam sequences:
- Rfam sequences: ~5 GB
- Total: ~8 GB
```

## Troubleshooting

### Rfam Download Fails

```bash
# Try alternative mirror
wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz

# Or download manually from:
# https://rfam.org/
```

### cmpress Fails

```bash
# Make sure Infernal is installed
which cmpress

# Check Rfam.cm file is valid
file Rfam.cm

# Should say: "data" or similar
```

### Out of Disk Space

```bash
# Check available space
df -h

# Use smaller database (Rfam sequences instead of RNAcentral)
# Or use external drive
```

### Slow Database Search

```bash
# Increase E-value threshold (less stringent)
# Edit pipeline code: evalue_threshold=0.1 instead of 0.01

# Or use fewer CPUs to reduce memory
--n_cpu 4
```

## Alternative: Use Pre-computed MSAs

If you can't download databases, you can:

1. **Use public MSAs**: Download pre-computed alignments from Rfam
2. **Request from collaborators**: Ask for MSAs for your specific RNAs
3. **Use single-sequence**: Pipeline works without databases (lower quality)

## Recommendations

### For Testing (Minimal Setup)
- Install external tools only
- No databases needed
- Use single-sequence fallback
- Good for verifying pipeline works

### For Production (Recommended Setup)
- Install external tools
- Download Rfam database (~2 GB)
- Optionally download sequence database (~5-50 GB)
- Best quality MSAs

### For HPC/Cluster
- Install tools in your environment
- Download databases to shared storage
- Use `--n_cpu` to match your allocation
- Consider batch processing

## Summary

**Minimum to run**: Nothing (uses fallback)
**Recommended**: Rfam database + external tools
**Optimal**: Rfam + sequence database + external tools

**Start with**: External tools + Rfam database
**Add later**: Sequence database if needed

See `QUICKSTART_RNA_MSA.md` for usage examples.
