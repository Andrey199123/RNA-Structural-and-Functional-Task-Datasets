# Task 5: RNA-RNA Interaction (RRI) Prediction

## Dataset
- **Name**: DeepMirTar
- **Source**: Wen et al., Bioinformatics 2018
- **Pairs**: 27,720 total (13,860 positive, 13,860 negative)
- **Status**: ✅ Downloaded

## Data Location
- `data/DeepMirTar_SdA/data/` - Pickle files with interaction data

## Dataset Details
- **miRNA sequences**: <26 nt
- **mRNA site sequences**: <53 nt (3' UTR regions)
- **Task**: Binary classification (interaction vs no interaction)
- **Balance**: 50% positive, 50% negative
- **Seed types**: Canonical and non-canonical

## Implementation Status
- [ ] Data preprocessing
- [ ] Extract train/val/test splits (CRITICAL - must match RNAErnie paper)
- [ ] Paired sequence handling (miRNA + mRNA)
- [ ] Model evaluation
- [ ] Results reporting
- [ ] Seed type analysis

## Critical Requirements
- **MUST use exact same train/val/test split as RNAErnie paper**
- This ensures fair comparison with baseline methods
- Handle paired sequences correctly (miRNA-mRNA pairs)

## Baseline
- DeepMirTar: AUROC ~0.90

## Notes
- Binary classification task
- Evaluation metrics: AUROC, AUPRC, Accuracy
- Seed type (canonical vs non-canonical) affects interaction strength
- miRNA-mRNA interaction is crucial for gene regulation
