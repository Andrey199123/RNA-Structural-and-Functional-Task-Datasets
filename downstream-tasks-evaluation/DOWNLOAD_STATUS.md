# Dataset Download Status

Last updated: May 24, 2026

## Summary

| Task | Dataset | Status | Location | Notes |
|------|---------|--------|----------|-------|
| Task 1 | RNAsnap2 (RSA) | ✅ Downloaded | `data/raw/task1_rsa/RNAsnap2/` | TR89, CN48, TS44 splits available |
| Task 2 | RNALocate v2.0 | ⏳ Manual Required | `data/raw/task2_subcellular/` | Under investigation - sequences may be too long |
| Task 3 | RNAcontact | ⏳ Manual Required | `data/raw/task3_contact_map/` | Need paper supplementary materials |
| Task 4 | PrismNet (RBP) | ✅ Downloaded | `data/raw/task4_rbp/PrismNet/data/clip_data/` | All 17 HeLa RBPs downloaded! |
| Task 5 | DeepMirTar (RRI) | ✅ Downloaded | `data/raw/task5_rri/DeepMirTar_SdA/data/` | miRNA-mRNA interaction data available |
| Task 6 | Secondary Structure | 🚫 On Hold | N/A | Awaiting approval - do not download |
| Task 7 | Long-Range BP | ⏳ Manual Required | `data/raw/task7_long_range/` | Need PDB Benchmark + bpRNA-1m TR0 |

## Detailed Status

### ✅ Task 1: RNA Solvent Accessibility (READY)

**Downloaded**: RNAsnap2 repository from GitHub

**Location**: `data/raw/task1_rsa/RNAsnap2/`

**Contents**:
- `benchmarking/sequences/`: 76 RNA sequence files
- `benchmarking/true_asa/`: Ground truth SASA values
- Test set IDs: `ts31_ids`, `ts45_ids`

**Next Steps**:
1. Identify TR89 training split (89 RNA chains)
2. Identify CN48 test split (48 protein-free RNAs)
3. Identify TS44 test split (44 protein-bound RNAs)
4. Parse FASTA sequences with position-wise labels
5. Implement per-nucleotide normalization

**Ready to implement**: YES ✅

---

### ✅ Task 4: RNA-Protein Binding (DOWNLOADED!)

**Downloaded**: PrismNet repository + full CLIP dataset from GitHub

**Location**: `data/raw/task4_rbp/PrismNet/data/clip_data/`

**Contents**:
- All 17 HeLa RBP datasets (TSV format):
  1. CSTF2_Hela.tsv
  2. EIF4A3_Hela.tsv
  3. ELAVL1_Hela.tsv
  4. HNRNPC_Hela.tsv
  5. HNRNPU_Hela.tsv
  6. METTL14_Hela.tsv
  7. METTL3_Hela.tsv
  8. PTBP1_Hela.tsv
  9. PTBP1PTBP2_Hela.tsv
  10. TIA1_Hela.tsv
  11. TIAL1_Hela.tsv
  12. U2AF2_Hela.tsv
  13. U2AF65_Hela.tsv
  14. UPF1_Hela.tsv
  15. WTAP_Hela.tsv
  16. YTHDF2_Hela.tsv
  17. eIF4AIII_Hela.tsv

**Download Command Used**:
```bash
cd data/raw/task4_rbp/PrismNet/data
wget https://zhanglabnet.oss-cn-beijing.aliyuncs.com/prismnet/data/clip_data.tgz
tar zxvf clip_data.tgz
```

**Next Steps**:
1. Parse TSV files to extract 101 nt sequences
2. Extract per-nucleotide binding labels
3. Implement train/val/test splits matching ERNIE-RNA paper
4. Optional: Extract icSHAPE features

**Ready to implement**: YES ✅

---

### ✅ Task 5: RNA-RNA Interaction (DOWNLOADED!)

**Downloaded**: DeepMirTar repository from GitHub

**Location**: `data/raw/task5_rri/DeepMirTar_SdA/data/`

**Contents**:
- `Mark_features_dic_all.pkl`: Feature dictionary for all miRNA-mRNA pairs
- `Mark_N_features_dic_all.pkl`: Negative features dictionary

**Repository**: https://github.com/Bjoux2/DeepMirTar_SdA

**Next Steps**:
1. Load pickle files to extract miRNA-mRNA pairs
2. Verify 27,720 total pairs (13,860 positive + 13,860 negative)
3. Extract sequences (miRNA <26 nt, mRNA <53 nt)
4. Implement train/val/test splits matching RNAErnie paper
5. Identify canonical vs non-canonical seed types

**Ready to implement**: YES ✅

---

### ⏳ Task 2: RNA Subcellular Localization (MANUAL REQUIRED)

**Status**: Under investigation - sequences may exceed model capacity (119-12,000 nt)

**Source**: http://www.rna-society.org/rnalocate/

**Manual Steps**:
1. Visit RNALocate website
2. Download RNALocate v2.0 dataset
3. Place in `data/raw/task2_subcellular/`
4. Analyze sequence length distribution
5. Identify reference papers (e.g., LOCAS)
6. Extract train/test splits

**Investigation Required**:
- Sequence length analysis
- Model capacity assessment
- Decision on task inclusion

---

### ⏳ Task 3: RNA Contact Map Prediction (MANUAL REQUIRED)

**Source**: Sun et al., Bioinformatics 2021 - RNAcontact dataset

**Expected**: 301 filtered sequences (TR221, TR168, VL43, TS80 splits)

**Manual Steps**:
1. Access paper: Sun et al., Bioinformatics 2021
2. Download supplementary materials
3. Look for RNAcontact dataset
4. Verify splits: TR221 (train), TR168+VL43 (val), TS80 (test)
5. Place in `data/raw/task3_contact_map/`

**CRITICAL**: Must match ERNIE-RNA paper splits exactly

---

### ⏳ Task 5: RNA-RNA Interaction (MANUAL REQUIRED)

**Source**: Wen et al., Bioinformatics 2018 - DeepMirTar dataset

**Expected**: 27,720 miRNA-mRNA pairs (13,860 positive + 13,860 negative)

**Manual Steps**:
1. Access paper: Wen et al., Bioinformatics 2018
2. Check supplementary materials or GitHub
3. Download DeepMirTar dataset
4. Verify balanced dataset (50% positive, 50% negative)
5. Place in `data/raw/task5_rri/`

**CRITICAL**: Must match RNAErnie paper train/val/test splits exactly

---

### 🚫 Task 6: RNA Secondary Structure (ON HOLD)

**Status**: DO NOT DOWNLOAD - awaiting approval

**Reason**: Model was already trained on secondary structure data

**Action**: Wait for explicit approval from Azam before proceeding

---

### ⏳ Task 7: Long-Range Base Pair Prediction (MANUAL REQUIRED)

**Source**: RNA-km paper (Section 3.4)

**Expected**:
- Dataset A: PDB Benchmark (532 sequences, 158 with pseudoknots)
- Dataset B: bpRNA-1m TR0 (10,814 sequences, 999 with pseudoknots)

**Manual Steps**:
1. Access RNA-km paper (Section 3.4)
2. Download PDB Benchmark dataset
3. Download bpRNA-1m TR0 split
4. Place in `data/raw/task7_long_range/`

**Requirements**:
- DSSR tool for structure annotation
- CD-HIT-EST for redundancy filtering (80% identity)

---

## Next Actions

### Immediate (Can Start Now)

1. **Task 1 (RSA)**: 
   - Inspect RNAsnap2 data format
   - Identify TR89, CN48, TS44 splits
   - Start implementing preprocessing pipeline

2. **Task 4 (RBP)**:
   - Parse TSV files for all 17 HeLa RBPs
   - Extract 101 nt sequences and binding labels
   - Implement preprocessing pipeline

3. **Task 5 (RRI)**:
   - Load DeepMirTar pickle files
   - Extract miRNA-mRNA pairs
   - Implement preprocessing pipeline

### High Priority (Manual Downloads Needed)

2. **Task 4 (RBP)**:
   - Check PrismNet README for full dataset download
   - Obtain all 17 RBP datasets

3. **Task 3 (Contact Map)**:
   - Download RNAcontact dataset from paper

4. **Task 5 (RRI)**:
   - Download DeepMirTar dataset from paper

5. **Task 7 (Long-Range)**:
   - Download PDB Benchmark and bpRNA-1m TR0

### Investigation Required

6. **Task 2 (Subcellular)**:
   - Download RNALocate v2.0
   - Analyze sequence lengths
   - Make inclusion decision

### Blocked

7. **Task 6 (Secondary Structure)**:
   - Awaiting approval from Azam
   - Do not proceed

---

## Implementation Priority

Based on data availability:

1. **Start with Task 1** - Data ready, can begin implementation
2. **Complete manual downloads** - Tasks 3, 4, 5, 7
3. **Investigate Task 2** - Analyze sequence lengths
4. **Wait on Task 6** - Blocked pending approval

---

## Tools Required

Before implementing preprocessing pipelines, ensure these tools are installed:

- **DSSR** (Task 7): RNA structure annotation
  - Download from: http://x3dna.org/
  
- **CD-HIT-EST** (Tasks 3, 7): Redundancy removal
  - Install: `conda install -c bioconda cd-hit`
  
- **BLASTclust** (Task 3): Alternative redundancy removal
  - Install: `conda install -c bioconda blast`
