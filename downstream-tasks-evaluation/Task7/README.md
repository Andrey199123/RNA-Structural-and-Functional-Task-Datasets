# Task 7: Long-Range Base Pair Prediction

## Datasets
This task uses **TWO separate datasets**:

### Dataset A: PDB Benchmark
- **Source**: 532 non-redundant RNA sequences from PDB
- **Processing**: 
  - Remove protein and DNA binding partners
  - CD-HIT-EST redundancy filtering at 80% identity
  - DSSR tool for secondary structure annotation
- **Pseudoknots**: 158 sequences contain pseudoknots
- **Status**: ⏳ Need to download

### Dataset B: bpRNA-1m TR0 Split
- **Source**: bpRNA-1m database (TR0 split in the original paper)
- **Sequences**: 10,814 sequences (NOT the full bpRNA-1m dataset)
- **Pseudoknots**: 999 sequences contain pseudoknots
- **Usage**: Zero-shot evaluation ONLY (no fine-tuning)
- **Status**: ✅ Downloaded via `multimolecule/bprna-spot-0` (TR0/VL0/TS0)

## Data Location
- `data/task7_long_range/` - Both datasets will be saved here

## How to Obtain Datasets

### Option 1: Hugging Face (Recommended for TR0)
```python
from datasets import load_dataset
ds = load_dataset("multimolecule/bprna-spot-0")
# ds["train"]      -> TR0 (10814)
# ds["validation"] -> VL0 (1300)
# ds["test"]       -> TS0 (1305)
```

### Option 2: bpRNA Website
- Visit: https://bprna.cgrb.oregonstate.edu/
- Look for TR0 split download (NOT the full 102,318 structures!)
- **Important**: Only download TR0 (10,814 seqs), not full bpRNA-1m (53.3 GB)

### Option 3: Contact Authors
- For PDB Benchmark (532 sequences)
- For TR0 split if not available online

## Dataset Details
- **Long-range definition**: |i-j| ≥ 101 (sequence separation)
- **Task**: Binary classification (base pair vs no base pair)
- **Focus**: Pseudoknots and long-range interactions
- **Evaluation**: Precision, Recall, F1-score

## Implementation Status
- [ ] Download PDB Benchmark (532 sequences)
- [x] Download bpRNA TR0/VL0/TS0 split (10,814/1,300/1,305)
- [ ] Apply CD-HIT-EST redundancy filtering (80% identity)
- [ ] Run DSSR for structure annotation
- [ ] Filter for long-range pairs (|i-j| ≥ 101)
- [ ] Separate pseudoknot-containing sequences
- [ ] Model evaluation (zero-shot on TR0)
- [ ] Results reporting

## Critical Requirements
- **TR0 is for zero-shot evaluation ONLY** - no fine-tuning allowed
- **Long-range pairs only**: |i-j| ≥ 101
- **External tools required**: CD-HIT-EST, DSSR

## Baseline
- Various methods: F1-score ~0.30-0.50 for long-range pairs
- Pseudoknots are particularly challenging

## Notes
- Long-range base pairs are harder to predict than local pairs
- Pseudoknots involve non-nested base pairing
- DSSR tool extracts base pairs from 3D structures
- CD-HIT-EST removes redundant sequences at 80% identity threshold
- **DO NOT download the full 53.3 GB bpRNA-1m dataset** - only need TR0 split!
- This repo stores the split files under `data/bprna_spot0/` (TR0/VL0/TS0 JSONL).
