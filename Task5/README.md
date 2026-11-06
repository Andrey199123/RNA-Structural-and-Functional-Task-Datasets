# miRNA–target interaction

Source: MirTarRAW/DeepMirTar as distributed by [RNAErnie](https://github.com/CatIIIIIIII/RNAErnie). Run `python Task5/prepare.py /path/to/Task5/splits` on the original project's CSVs. Each adjacent miRNA and target record shares a pair ID and label.

The 27,719 pairs retain the local 80/10/10 partition and its 26/53-nt padded inputs. The local validation and test partitions together correspond to the held-out portion, so the test partition alone is not the full RNAErnie test set.
