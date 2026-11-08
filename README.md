# RNA Structural and Functional Task Datasets

Curated preparation scripts and reproducible dataset snapshots for the RNAtomy
RNA structure and function study.

This repository turns the team's final local dataset exports into small,
reviewable task folders. Each folder contains the preparation code, the
resulting sequence files, and the companion labels or structural annotations
needed to interpret them. It contains **dataset preparation only**. It does not
contain model training, prediction code, checkpoints, or reported results.

## Dataset overview

| Folder | Application | Included benchmark material |
|:--|:--|:--|
| [`Task1`](Task1) | Solvent accessibility | TR95 / VL24 / TS44, 95 / 24 / 44 RNAs |
| [`Task2`](Task2) | IRES identification | 12,983 training / 1,164 test sequences |
| [`Task3`](Task3) | 3D contact maps | 168 training / 43 validation / 80 test RNAs |
| [`Task4`](Task4) | RNA-protein binding | 17 HeLa RBP datasets, 101-nt windows |
| [`Task5`](Task5) | miRNA-target interaction | 27,719 labeled pairs |
| [`Task7`](Task7) | Long-range base pairs | TR0 / VL0 / TS0, 10,814 / 1,300 / 1,305 sequences |

The folder numbers follow the applications in the current paper draft. Task 6,
mean ribosome load, is intentionally omitted because the exact final experiment
splits were not available for export. The folders are not the old checkout's
localization or inverse-folding projects.

## Quick start

Requirements: Python 3.11 or newer.

```bash
python -m pip install -r requirements.txt
python Task1/prepare.py --output /tmp/task1
```

Replace `Task1` and `/tmp/task1` with any populated task folder and output
directory. Every included task has a bundled `inputs.zip`, so the command works
without the old workspace, downloads, model weights, or machine-specific paths.

Each script also accepts an external input directory or ZIP archive when a
different source export needs to be processed. See the README inside each task
folder for its source, split details, and format-specific notes.

## What is preserved

FASTA files carry the sequences and identifiers. Companion files preserve the
information that FASTA cannot represent on its own:

- RSA labels for solvent accessibility
- contact and validity matrices for 3D contact prediction
- icSHAPE tracks for RBP binding
- dot-bracket structures for long-range base-pair prediction

The bundled archives are labeled local-input snapshots used to reproduce the
included exports. They are not presented as raw upstream releases or as a claim
that every local split is identical to the final paper experiment. The per-task
READMEs document these distinctions rather than silently changing the data.

## Reproducibility notes

- Preparation scripts are deterministic for the included snapshots.
- Scripts validate input structure and write outputs to a user-selected folder.
- Scripts do not download data, fold structures, or require access to FIU servers.
- The repository contains no temporary files, caches, checkpoints, or unrelated
  RNA-protein baseline materials.

## Contributors

Dataset preparation was a team effort under Azam Shirali's guidance. This list
acknowledges project participation and does not claim exclusive authorship of
every script or dataset.

- [Andrey Vasilyev](https://github.com/Andrey199123)
- [Ishaan Garg](https://github.com/ishaantek)
- [Pavitra Soni](https://github.com/pavitrasoni113)

## Sources and use

Every task README links to the relevant upstream dataset or paper. Please cite
the original dataset and method papers when using these files. Upstream data
remain subject to their authors' terms; this repository does not assert a blanket
license over third-party datasets. Check redistribution terms before publishing
mirrors or incorporating the files into another release.
