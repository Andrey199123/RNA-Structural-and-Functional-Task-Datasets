# Task 4: RNA-Protein Binding (RBP) Prediction

## Dataset
- **Name**: PrismNet (17 HeLa RBPs)
- **Source**: github.com/kuixu/PrismNet
- **RBPs**: 17 RNA-binding proteins in HeLa cell environment
- **Sequences**: 101 nt fixed-length (1,827 to 15,002 per RBP)
- **Status**: ✅ Downloaded

## Data Location
- `data/PrismNet/data/clip_data/` - TSV files for all 17 RBPs

## RBP List
1. CSTF2
2. EIF4A3 (eIF4AIII)
3. ELAVL1
4. HNRNPC
5. HNRNPU
6. METTL14
7. METTL3
8. PTBP1
9. PTBP1PTBP2
10. TIA1
11. TIAL1
12. U2AF2 (U2AF65)
13. UPF1
14. WTAP
15. YTHDF2

## Dataset Details
- **Format**: TSV files with sequences and binding labels
- **Sequence length**: Fixed 101 nt
- **Task**: Binary classification (binding vs non-binding)
- **Optional**: icSHAPE features available for comparison

## Implementation Status
- [ ] Data preprocessing
- [ ] Extract train/val/test splits (CRITICAL - must match ERNIE-RNA paper)
- [ ] Model evaluation per RBP
- [ ] Results reporting
- [ ] Optional: Compare with/without icSHAPE features

## Critical Requirements
- **MUST use exact same train/val/test split as ERNIE-RNA paper**
- This ensures fair comparison with baseline methods
- Each RBP is evaluated separately

## Baseline
- PrismNet: AUROC ~0.75-0.85 (varies by RBP)

## Notes
- Binary classification task for each RBP
- Evaluation metrics: AUROC, AUPRC
- icSHAPE features provide in vivo RNA structure information
