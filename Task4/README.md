# RNA–protein binding

Reproduce the included FASTAs and labels without the old project:

```bash
python Task4/prepare.py --output /tmp/task4
```

The default `inputs.zip` is the minimal labeled local-input snapshot used by
this export script. It preserves the existing sequence order and partitions.
The source-specific command below is optional, not required for reproduction.

Source: [PrismNet](https://github.com/kuixu/PrismNet), 17 HeLa datasets. Run `python Task4/prepare.py /path/to/Task4/hdf5` on the team's 72/8/20 HDF5 exports. The local split uses seed 42; matching proportions does not establish exact ERNIE-RNA membership.

FASTA headers contain the RBP dataset, zero-based row index, and binary label. The matching NPZ preserves per-base icSHAPE measurements by RBP and row.
