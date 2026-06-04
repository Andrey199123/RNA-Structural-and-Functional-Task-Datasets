# Task 3: RNA Contact Map Prediction

## Dataset
- **Name**: RNAcontact
- **Source**: Sun et al., Bioinformatics 2021 (doi: 10.1093/bioinformatics/btaa1005)
- **Splits**: TR221 (training), TR168 + VL43 (validation), TS80 (test)
- **Total**: 301 filtered sequences from 1,786 PDB RNA structures
- **Status**: ✅ Downloaded + preprocessed

## Data Location
- `data/` - Raw RNAcontact benchmark archives + extracted files
- `preprocessed/` - Per-sequence `.npz` files with contact maps

## How to Obtain Dataset
The public benchmark archives are available from the RNAcontact benchmark page.

## Dataset Details
- **Filtering**: 301 sequences selected from 1,786 PDB RNA structures
- **Contact definition**: Minimum atom-pair distance ≤ 8Å
- **Long-range contacts**: |i-j| ≥ 24 (sequence separation)
- **Format**: Sequences + contact matrices

## Implementation Status
- [x] Download and extract data
- [x] Contact map generation
- [x] Data preprocessing
- [ ] Model evaluation
- [ ] Results reporting

## Critical Requirements
- **MUST use exact splits from ERNIE-RNA paper** (TR221, TR168, VL43, TS80)
- This ensures fair comparison with baseline methods

## Baseline
- RNAcontact: Precision ~0.45-0.50 for long-range contacts

## Notes
- Contact maps are binary matrices indicating 3D proximity
- Long-range contacts (|i-j| ≥ 24) are more challenging to predict
- Evaluation metric: Precision at top L/5, L/2, L predictions (L = sequence length)
