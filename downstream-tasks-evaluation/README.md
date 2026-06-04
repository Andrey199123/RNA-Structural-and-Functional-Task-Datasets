# Downstream Tasks Evaluation

This folder contains task-specific scripts, downloaded datasets, preprocessing outputs, and run reports for the `downstream_tasks` orchestrator.

## Directory Structure

```
downstream-tasks-evaluation/
├── Task1/..Task7/              # Per-task data + preprocessed outputs
├── config/                     # task_configs.yaml, preprocessing_configs.yaml, validation_specs.yaml
├── logs/                       # task logs
└── reports/                    # per-task status.json + consolidated_status.json + reproducibility.json
```

## Quick Start

### 1) Create a venv (recommended on macOS)

If `pip` shows `error: externally-managed-environment`, create a virtual environment and install packages there:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -r ../downstream_tasks/requirements.txt
```

### 2) Run the orchestrator

```bash
python -m downstream_tasks.cli run-all
python -m downstream_tasks.cli status
python -m downstream_tasks.cli consolidate
```

Common options:

```bash
python -m downstream_tasks.cli run-all --parallel --max-workers 4
python -m downstream_tasks.cli run-task task7_long_range --force
```

## Task Notes

- **Task 1 (RSA)**: uses `Task1/data/rnasol_benchmark/` (TR89 / TS44 / CN48). `prepare_task1_rsa.py` extracts the local tarballs if needed; no network download.
- **Task 3 (Contact Map)**: downloads the public RNAcontact benchmark archives and computes 8Å long-range contacts (|i−j| ≥ 24).
- **Task 7 (Long-Range)**: uses `multimolecule/bprna-spot-0` TR0 exported by `Task7/download_bprna_tr0.py`.
- **Task 6 (Secondary Structure)**: disabled by default (see `config/task_configs.yaml`) pending explicit approval.

## Kiro Specs

The original spec documents live in `../.kiro/specs/downstream-tasks-evaluation/`.
