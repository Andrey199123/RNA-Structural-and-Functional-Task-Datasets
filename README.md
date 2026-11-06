# RNA Structural and Functional Task Datasets

Preparation scripts and sequence datasets assembled for RNAtomy. This repository
contains dataset preparation only, not model training, prediction, or results.

| Folder | Application | Included sequences |
| --- | --- | --- |
| Task1 | Solvent accessibility | TR95 / VL24 / TS44: 95 / 24 / 44 |
| Task2 | IRES identification | 12,983 train / 1,164 test |
| Task3 | 3D contact maps | 168 train / 43 validation / 80 test |
| Task4 | RBP binding | 17 HeLa datasets; 101-nt windows |
| Task5 | miRNA–target interaction | 27,719 pairs |
| Task7 | Long-range base pairs | TR0 / VL0 / TS0: 10,814 / 1,300 / 1,305 |

**This is a review repository, not yet an assertion that every local export is
identical to the final paper's experiment files.** Mean ribosome load is excluded
because the exact final experiment splits are unavailable.
The per-task READMEs identify the upstream sources and remaining split questions.

## Using the files

Each task has one `prepare.py`, FASTA files where source data are available,
and a short source README. Each populated task includes an `inputs.zip` snapshot
of the labeled inputs needed to reproduce its published files. These are the
team's local exports, not a claim of raw upstream data or final-paper approval.
They do not require the old workspace, old pipeline code, or
models. Small companion files retain labels and structural measurements that
FASTA cannot represent: RSA tables, contact maps, icSHAPE tracks, or dot-bracket
annotations. Never discard those files when reproducing a labeled benchmark.

Install Python 3.11 or newer and `python -m pip install -r requirements.txt`.
For Tasks 1–5 and 7, run `python Task1/prepare.py --output /tmp/task1` (substitute
the task number) to reproduce from the bundled input archive, without downloads.
An optional input-directory or ZIP argument still supports external source files.
No MRL data or placeholder preparation script is included. No machine-specific paths are embedded.
Scripts never download or fold structures.

## Folder numbering

The folder numbers follow the applications in the current paper draft.
Task2 is IRES; Task6 (MRL) is omitted. They are **not** the old checkout's Task2
localization and Task6 inverse-folding projects; both old projects are excluded.

## Contributors and upstream data

- [Andrey Vasilyev](https://github.com/Andrey199123)
- [Ishaan Garg](https://github.com/ishaantek)
- [Pavitra Soni](https://github.com/pavitrasoni113)

Dataset preparation was a team effort under Azam Shirali's guidance. This list
acknowledges participation, not exclusive authorship of every script or dataset.
The new export scripts organize the team's existing local outputs. This is the
original research repository, renamed RNA-Structural-and-Functional-Task-Datasets.

Upstream datasets remain subject to their authors' terms and citations; no
blanket license over third-party data is asserted here. Review redistribution
terms and confirm final experiment memberships before making the repository
public or linking it as the paper's final dataset release. The original private
research repository is preserved separately; no RNA–protein baseline experiments,
alignments, model checkpoints, caches, temporary files, or unrelated papers
are included here.
