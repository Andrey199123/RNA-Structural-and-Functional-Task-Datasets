# RNA 3D contact maps

Reproduce the included FASTAs and labels without the old project:

```bash
python Task3/prepare.py --output /tmp/task3
```

The default `inputs.zip` is the minimal labeled local-input snapshot used by
this export script. It preserves the existing sequence order and partitions.
The source-specific command below is optional, not required for reproduction.

Source: [RNAcontact benchmark](https://yanglab.qd.sdu.edu.cn/RNAcontact/benchmark/). Run `python Task3/prepare.py /path/to/baseline_comparisons/work/contact_c1` on the original project's coordinate-derived JSONL export (168/43/80 RNAs).

The local export uses C1′ distances strictly below 8 Å, not legacy heavy-atom labels. Each FASTA `matrix` index identifies `contact_N` and `valid_N` in the matching NPZ. Missing coordinates are masked. Confirm these local partitions and masks match the final RNAtomy run.
