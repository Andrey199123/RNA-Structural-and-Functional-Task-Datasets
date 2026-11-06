# Mean ribosome load

Source: Sample et al. 5′ UTR assays ([GEO GSE114002](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114002)); the paper draft names the Random7600 and Human7600 evaluations.

**Pending the exact experiment files.** The original checkout does not contain verified training or 7,600-example evaluation memberships. Do not resample GEO and call it the same benchmark. Once the exact labeled CSVs are available, name them `train.csv`, `Random7600.csv`, and `Human7600.csv`, with `sequence,mrl` columns, and run `python Task6/prepare.py /path/to/exact_splits`.
