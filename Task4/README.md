# RNA–protein binding

Source: [PrismNet](https://github.com/kuixu/PrismNet), 17 HeLa datasets. Run `python Task4/prepare.py /path/to/Task4/hdf5` on the team's 72/8/20 HDF5 exports. The local split uses seed 42; matching proportions does not establish exact ERNIE-RNA membership.

FASTA headers contain the RBP dataset, zero-based row index, and binary label. The matching NPZ preserves per-base icSHAPE measurements by RBP and row.
