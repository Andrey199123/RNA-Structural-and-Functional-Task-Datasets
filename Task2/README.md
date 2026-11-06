# Internal ribosome entry site identification

Reproduce the included FASTAs and labels without the old project:

```bash
python Task2/prepare.py --output /tmp/task2
```

The default `inputs.zip` is the minimal labeled local-input snapshot used by
this export script. It preserves the existing sequence order and partitions.
The source-specific command below is optional, not required for reproduction.

Source: [structRFM task-data release](https://zenodo.org/records/16754363), `IRES/data/ires_seq`. Run `python Task2/prepare.py /path/to/ires_seq` after extracting the upstream archive.

The FASTA headers carry the binary labels. The release contains 12,983 training and 1,164 test segments of 174 nt; the test set has 582 examples per class. No new validation split is imposed. This folder is IRES, not the postponed localization task from the old checkout.
