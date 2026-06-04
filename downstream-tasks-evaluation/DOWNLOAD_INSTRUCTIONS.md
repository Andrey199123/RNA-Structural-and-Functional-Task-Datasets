# Dataset Download Instructions

## Summary

This document provides instructions for downloading the remaining datasets for the 7 downstream tasks.

## ✅ Already Downloaded (3/7 tasks)

1. **Task 1 (RSA)**: RNAsnap2 - Downloaded ✓
2. **Task 4 (RBP)**: PrismNet (17 HeLa RBPs) - Downloaded ✓
3. **Task 5 (RRI)**: DeepMirTar - Downloaded ✓

## ⏳ Need to Download (3/7 tasks)

### Task 2: RNA Subcellular Localization (UNDER INVESTIGATION)

**Status**: Downloaded but under investigation (sequences may be too long: 119-12,000 nt)

**Preferred file location (manual download)**: `data/raw/task2_subcellular/RNALocate_v2.0.txt`

**Legacy file (annotations only; typically no sequences)**: `Task2/data/All RNA subcellular localization information.txt`

**Next steps**:
1. Analyze sequence length distribution
2. Determine if sequences exceed model capacity
3. Decide on task inclusion

### Task 3: RNA Contact Map Prediction

**Dataset**: RNAcontact (TR221, TR168, VL43, TS80 splits)

**How to obtain**:
1. Email the authors: zengjy321@tsinghua.edu.cn
2. Use the email template in `EMAIL_TEMPLATE_RNACONTACT.txt`
3. Request the TR221, TR168, VL43, and TS80 splits with sequences and contact maps

**Alternative**: The dataset is NOT publicly available on GitHub, Zenodo, or Figshare based on our search. Emailing the authors is the recommended approach.

### Task 7: Long-Range Base Pair Prediction

**Datasets needed**:
- **Dataset A**: PDB Benchmark (532 sequences)
- **Dataset B**: bpRNA-1m TR0 split (10,814 sequences)

**IMPORTANT**: Only download the TR0 split, NOT the full 53.3 GB bpRNA-1m dataset!

**How to obtain**:

#### Option 1: Python script (Recommended)
```bash
cd Task7
python3 download_bprna_tr0.py
```

This script downloads only the TR0 split using Hugging Face datasets library.

#### Option 2: Hugging Face manually
```python
from datasets import load_dataset
ds = load_dataset("multimolecule/bprna-spot-0")
# ds["train"]      -> TR0 (10,814)
# ds["validation"] -> VL0 (1,300)
# ds["test"]       -> TS0 (1,305)
```

#### Option 3: bpRNA website
- Visit: https://bprna.cgrb.oregonstate.edu/
- Look for TR0 split download
- **DO NOT download the full 102,318 structures (53.3 GB)!**

#### Option 4: Contact authors
- For PDB Benchmark (532 sequences)
- For TR0 split if not available through other methods

## 🚫 On Hold (1/7 tasks)

### Task 6: RNA Secondary Structure Prediction

**Status**: BLOCKED - Awaiting explicit approval from Azam

**Reason**: Model was trained on secondary structure data, so evaluation may not be meaningful.

**Action**: DO NOT DOWNLOAD until you receive explicit approval.

## Quick Reference

| Task | Dataset | Status | Action |
|------|---------|--------|--------|
| 1 | RNAsnap2 | ✅ Downloaded | Ready for implementation |
| 2 | RNALocate v2.0 | ⚠️ Under investigation | Analyze sequence lengths |
| 3 | RNAcontact | ⏳ Need to request | Email authors |
| 4 | PrismNet | ✅ Downloaded | Ready for implementation |
| 5 | DeepMirTar | ✅ Downloaded | Ready for implementation |
| 6 | Secondary Structure | 🚫 Blocked | Awaiting approval |
| 7 | bpRNA-1m TR0 + PDB | ⏳ Need to download | Run download script |

## Implementation Priority

Based on Azam's request to send completed tasks incrementally:

1. **Start immediately**: Tasks 1, 4, 5 (data ready)
2. **After download**: Task 7 (download TR0 split)
3. **After receiving data**: Task 3 (email authors)
4. **After investigation**: Task 2 (analyze sequence lengths)
5. **Blocked**: Task 6 (awaiting approval)

## Notes

- Each Task folder (Task1-Task7) contains a README.md with detailed information
- See `../. kiro/specs/downstream-tasks-evaluation/tasks.md` for the full implementation plan
- CRITICAL: Tasks 3, 4, 5 must use exact paper splits for fair comparison
- Task 7 TR0 split is for zero-shot evaluation only (no fine-tuning)
