# RNA solvent accessibility

Reproduce the included FASTAs and labels without the old project:

```bash
python Task1/prepare.py --output /tmp/task1
```

The default `inputs.zip` is the minimal labeled local-input snapshot used by
this export script. It preserves the existing sequence order and partitions.
The source-specific command below is optional, not required for reproduction.

[RNAsnap2](https://doi.org/10.1093/bioinformatics/btaa652) TR95 (train) and VL24 (validation), and the [RNAsol](https://doi.org/10.1093/bioinformatics/bty876) TS44 test labels from its TINKER file. The local export contains 95/24/44 sequences, respectively; TS45 is not included.

Run `python Task1/prepare.py /path/to/baseline_comparisons/work/rsa` on the original project's prepared JSONL directory. The compressed tables retain normalized accessibility at one-based nucleotide positions; omitted positions are masked. The local conversion uses maxima of 400 Å² for A/G/N and 350 Å² for C/U. Confirm the normalization against the final paper run before comparing scores.
