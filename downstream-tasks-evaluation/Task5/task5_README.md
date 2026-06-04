# Task 5 — RNA-RNA Interaction Prediction

## Source
DeepMirTar (Wen et al., Bioinformatics 2018)
Dataset: MirTarRAW (DeepMirTar + miRAW combined)
As used by: RNAErnie (Wang et al., Nature Machine Intelligence 2024)
  github.com/CatIIIIIIII/RNAErnie

## Task
Binary classification: predict whether a miRNA and an mRNA target site interact.

## Dataset stats
- Total pairs:               27,719
- Positive (interacting):    13,859
- Negative (non-interacting):13,860

## Sequence format
- miRNA (a_seq): padded/trimmed to 26 nt with N, T→U converted
- mRNA  (b_seq): padded/trimmed to 53 nt with N, T→U converted
- label: 0 = no interaction, 1 = interaction

## Split — RNAErnie convention (rna_rna_interaction.py)
- Shuffle: numpy.random.RandomState(seed=0)
- Train: first 80%  → 22,175 pairs  (identical to RNAErnie train)
- Val:   next  10%  → 2,772 pairs
- Test:  last  10%  → 2,772 pairs   (comparable to RNAErnie results)

## Files
  splits/train.csv
  splits/val.csv
  splits/test.csv

## CSV columns
  a_name  miRNA identifier
  a_seq   miRNA sequence (26 nt, N-padded, RNA)
  b_name  mRNA identifier
  b_seq   mRNA target site (53 nt, N-padded, RNA)
  label   0 or 1

## Load example
import pandas as pd
train = pd.read_csv("splits/train.csv")
X_mirna = train["a_seq"].values
X_mrna  = train["b_seq"].values
y       = train["label"].values
