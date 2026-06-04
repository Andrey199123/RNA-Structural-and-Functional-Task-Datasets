# Downstream Tasks Dataset Sources

This document lists all dataset sources for the 7 downstream tasks and provides download instructions.

## Task 1: RNA Solvent Accessibility (RSA) Prediction

**Status**: Ready for download

**Dataset**: RNAsnap2 benchmark (TR89, CN48, TS44)
- TR89: 89 RNA chains (training)
- CN48: 48 protein-free RNAs (testing)
- TS44: 44 protein-bound RNAs (testing)

**Source**: 
- Paper: Yang et al., RNA 2017; Sun et al., Bioinformatics 2019
- Repository: https://github.com/jaswindersingh2/RNAsnap2

**Download Method**: Git clone (automated in download_datasets.py)

**Notes**:
- Labels already provided (no FreeSASA computation needed)
- FASTA sequences with position-wise float labels
- Evaluation metric: Pearson correlation

---

## Task 2: RNA Subcellular Localization

**Status**: UNDER INVESTIGATION - sequences may be too long (119-12,000 nt)

**Dataset**: RNALocate v2.0 (~3,300 entries after deduplication)

**Source**:
- Paper: Cui et al., Nucleic Acids Research 2022
- Website: http://www.rna-society.org/rnalocate/
- Download: rnalocate.org and supplementary NAR materials

**Download Method**: Manual download required

**Manual Steps**:
1. Visit http://www.rna-society.org/rnalocate/
2. Navigate to download section
3. Download RNALocate v2.0 dataset
4. Place files in `data/raw/task2_subcellular/`

**Investigation Required**:
- Analyze sequence length distribution
- Identify papers using this dataset (e.g., LOCAS - Rahman Abir et al., 2025)
- Extract train/test splits from reference papers
- Determine if sequences can be truncated or task should be excluded

---

## Task 3: RNA Contact Map Prediction

**Status**: Ready for download (manual)

**Dataset**: RNAcontact (301 filtered sequences)
- Initial: 1,786 RNA sequences from PDB
- Final: 301 sequences (>5 contacts)
- Splits: TR221 (training), TR168 + VL43 (validation), TS80 (testing)

**Source**:
- Paper: Sun et al., Bioinformatics 2021
- Dataset: Leontis and Zirbel PDB dataset

**Download Method**: Manual download from paper supplementary materials

**Manual Steps**:
1. Access paper: Sun et al., Bioinformatics 2021
2. Download supplementary materials or contact authors
3. Look for RNAcontact dataset with TR221/TR168/VL43/TS80 splits
4. Place files in `data/raw/task3_contact_map/`

**CRITICAL**: Must use exact splits from ERNIE-RNA paper for fair comparison

**Notes**:
- Contact definition: minimum atom-pair distance < 8Å
- Evaluate only long-range contacts (|i-j| ≥ 24)
- Evaluation: Long-Range Top-L precision (L/10, L/5, L/2, L/1)

---

## Task 4: RNA-Protein Binding Prediction

**Status**: Ready for download

**Dataset**: PrismNet dataset (17 RBPs in HeLa cells)
- 17 RBPs in HeLa environment
- Fixed sequence length: 101 nt
- Sequences per RBP: 1,827 to 15,002

**Source**:
- Paper: Xu et al., Nucleic Acids Research 2023
- Repository: https://github.com/kuixu/PrismNet

**Download Method**: Git clone (automated in download_datasets.py)

**CRITICAL**: Must use exact same train/val/test split as ERNIE-RNA paper

**Notes**:
- HDF5 (.h5) file format
- Optional icSHAPE features for comparison
- Evaluation: Mean AUPRC across all 17 proteins

---

## Task 5: RNA-RNA Interaction

**Status**: Ready for download (manual)

**Dataset**: DeepMirTar (27,720 pairs)
- 13,860 positive pairs
- 13,860 negative pairs
- miRNA: <26 nt
- mRNA: <53 nt (3' UTR)

**Source**:
- Paper: Wen et al., Bioinformatics 2018
- Dataset: DeepMirTar

**Download Method**: Manual download (contact authors or check paper supplementary)

**Manual Steps**:
1. Access paper: Wen et al., Bioinformatics 2018
2. Check supplementary materials or GitHub
3. Download DeepMirTar dataset
4. Place files in `data/raw/task5_rri/`

**CRITICAL**: Must use exact same train/val/test split as RNAErnie paper

**Notes**:
- Tab-separated interaction records
- Balanced dataset (50% positive, 50% negative)
- Evaluation: Accuracy, Recall, Precision, F1, AUC

---

## Task 6: RNA Secondary Structure Prediction

**Status**: ON HOLD - DO NOT DOWNLOAD

**Reason**: Model was already trained on secondary structure data - evaluation may not be meaningful

**Datasets** (if approved):
- bpRNA-1m (TR0/TV0/TS0 splits)
  - TR0: 10,814 training sequences
  - VL0: 1,300 validation sequences
  - TS0: 1,305 testing sequences
- RNAStralign (37,149 structures)
- ArchiveII (3,975 structures)
- Szikszai 2022 Cross-Family Split

**DO NOT PROCEED** until explicit approval from Azam

---

## Task 7: Long-Range Base Pair Prediction

**Status**: Ready for download (manual)

**Datasets**:

**Dataset A: PDB Benchmark**
- 532 non-redundant RNA sequences
- 158 sequences with pseudoknots
- Redundancy filtering: CD-HIT-EST at 80% identity

**Dataset B: bpRNA-1m TR0**
- 10,814 sequences
- 999 sequences with pseudoknots
- Zero-shot evaluation only (no fine-tuning)

**Source**:
- Paper: RNA-km paper (Section 3.4)
- Datasets: PDB + bpRNA-1m

**Download Method**: Manual download from paper supplementary materials

**Manual Steps**:
1. Access RNA-km paper (Section 3.4)
2. Download PDB Benchmark dataset (532 sequences)
3. Download bpRNA-1m TR0 split (10,814 sequences)
4. Place files in `data/raw/task7_long_range/`

**Notes**:
- Filter base pairs: |i-j| ≥ 101 for long-range evaluation
- Requires DSSR tool for structure annotation
- Evaluation: Accuracy for long-range pairs

---

## Download Priority

Based on Azam's request to send completed tasks incrementally:

1. **Task 1 (RSA)** - Automated download available
2. **Task 4 (RBP)** - Automated download available
3. **Task 3 (Contact Map)** - Manual download required
4. **Task 5 (RRI)** - Manual download required
5. **Task 7 (Long-Range)** - Manual download required
6. **Task 2 (Subcellular)** - Under investigation
7. **Task 6 (Secondary Structure)** - On hold

---

## Automated vs Manual Downloads

**Automated** (via download_datasets.py):
- Task 1: RNAsnap2 (GitHub)
- Task 4: PrismNet (GitHub)

**Manual Required**:
- Task 2: RNALocate v2.0 (website)
- Task 3: RNAcontact (paper supplementary)
- Task 5: DeepMirTar (paper supplementary)
- Task 7: PDB Benchmark + bpRNA-1m TR0 (paper supplementary)

**On Hold**:
- Task 6: All datasets (awaiting approval)

---

## External Tools Required

Some tasks require external bioinformatics tools:

1. **DSSR** (Task 7): RNA structure annotation
   - Download from: http://x3dna.org/
   - Version: ≥2.0

2. **CD-HIT-EST** (Task 3, Task 7): Redundancy removal
   - Install: `conda install -c bioconda cd-hit`
   - Version: ≥4.8

3. **BLASTclust** (Task 3): Alternative redundancy removal
   - Install: `conda install -c bioconda blast`
   - Version: ≥2.10

---

## Next Steps

1. Run `python download_datasets.py` to download automated datasets
2. Manually download remaining datasets following instructions above
3. Verify all datasets are in correct directories
4. Inspect data formats before implementing preprocessing pipelines
5. For Task 2: Perform investigation on sequence lengths before proceeding
