# Task 4 — RNA-Protein Binding Prediction

## Source
PrismNet (Xu et al., Nucleic Acids Research 2023)
https://github.com/kuixu/PrismNet

## Split (ERNIE-RNA convention)
- Train: 80%
- Val:   10%
- Test:  10%
- Random seed: 42

## Cell line
HeLa only

## HDF5 file structure
Each file hdf5/<RBP>_Hela.h5 has three groups: train / val / test
Each group contains:
  sequences  (N,)         raw 101-nt sequences as bytes
  one_hot    (N, 101, 5)  one-hot encoded ACGUN
  icshape    (N, 101)     icSHAPE scores, -1 = unmeasured
  labels     (N,)         0 = non-binding, 1 = binding

## Load example
import h5py
with h5py.File("hdf5/TIA1_Hela.h5", "r") as f:
    X = f["test/one_hot"][:]   # shape (N, 101, 5)
    y = f["test/labels"][:]    # shape (N,)

## RBPs processed (17)
HNRNPC, HNRNPU, CSTF2, METTL14, METTL3, TIA1, TIAL1, U2AF2, U2AF65, UPF1, WTAP, YTHDF2, EIF4A3, ELAVL1, PTBP1, PTBP1PTBP2, eIF4AIII
