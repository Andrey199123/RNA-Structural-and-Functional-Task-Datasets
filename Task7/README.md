# Long-range base-pair prediction

Reproduce the included FASTAs and labels without the old project:

```bash
python Task7/prepare.py --output /tmp/task7
```

The default `inputs.zip` is the minimal labeled local-input snapshot used by
this export script. It preserves the existing sequence order and partitions.
The source-specific command below is optional, not required for reproduction.

Source: [SPOT-RNA](https://github.com/jaswindersingh2/SPOT-RNA) bpRNA-1m partitions, via [bpRNA-spot-0](https://huggingface.co/datasets/multimolecule/bprna-spot-0). Run `python Task7/prepare.py /path/to/bprna_spot0` on the original project's JSONL exports.

TR0/VL0/TS0 contain 10,814/1,300/1,305 sequences. Dot-bracket files preserve the source pseudoknot symbols. Split names do not imply a training protocol: the draft also reports a TR0 evaluation, which must not be confused with training on TR0 and testing on TS0. No 532-structure substitute is included.
