#!/usr/bin/env python3
"""
Task 2: mRNA Subcellular Localization (RNALocate / DM3Loc benchmark) — data preparation.

This script prepares the DM3Loc/LOCAS-style multi-label mRNA localization benchmark:
  - Labels embedded in FASTA headers as a 7-bit one-hot string
  - We predict 6 compartments (exclude cytoplasm as an auxiliary class used in DM3Loc training)
  - Apply DM3Loc-style truncation to a fixed maximum length of 8000 nt:
      if len(seq) > 8000: seq[:4000] + seq[-4000:]
      else: keep whole sequence (no padding stored; pad during batching as needed)
  - Write gzipped JSONL records and deterministic 5-fold CV splits.

Input FASTA (default):
  Task2/data/dm3loc/modified_multilabel_seq_nonredundent.fasta

Outputs:
  Task2/preprocessed/dm3loc_trunc8000.jsonl.gz
  Task2/preprocessed/folds.json
  Task2/preprocessed/summary.json
  Task2/preprocessed/README.md
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

import numpy as np


THIS_DIR = Path(__file__).resolve().parent
DEFAULT_FASTA = THIS_DIR / "data" / "dm3loc" / "modified_multilabel_seq_nonredundent.fasta"
OUTPUT_DIR = THIS_DIR / "preprocessed"

LABEL_NAMES_6 = ["Nucleus", "Exosome", "Cytosol", "Ribosome", "Membrane", "ER"]
# DM3Loc 7-bit label order (from DM3Loc code / class weights):
#   0 Nucleus
#   1 Exosome
#   2 Cytosol
#   3 Cytoplasm (auxiliary; skipped in prediction output)
#   4 Ribosome
#   5 Membrane
#   6 Endoplasmic reticulum (ER)
IDX_6 = [0, 1, 2, 4, 5, 6]
IDX_CYTOPLASM = 3


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: str | None = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def parse_dm3loc_header(header: str) -> Dict[str, str | List[int]]:
    """
    Header format:
      0100000,ACCNUM:NM_001507,Gene_ID:2862,Gene_Name:MLNR
    """
    parts = header.split(",")
    bits = parts[0].strip()
    if not re.fullmatch(r"[01]{7}", bits):
        raise ValueError(f"Unexpected label bits: {bits!r}")
    meta: Dict[str, str | List[int]] = {"label_bits_7": [int(c) for c in bits]}
    for p in parts[1:]:
        if ":" not in p:
            continue
        k, v = p.split(":", 1)
        meta[k.strip()] = v.strip()
    meta["header"] = header
    return meta


def dm3loc_truncate(seq: str, *, max_len: int, left: int, right: int) -> Tuple[str, str]:
    seq = str(seq).upper().replace("T", "U")
    if len(seq) <= max_len:
        return seq, "none"
    if left + right != max_len:
        raise ValueError("left+right must equal max_len for end-concat truncation")
    return (seq[:left] + seq[-right:]), f"end_concat_{left}_{right}"


def compute_length_summary(lengths: np.ndarray) -> Dict[str, float]:
    if lengths.size == 0:
        return {"n": 0.0, "min": 0.0, "p50": 0.0, "p90": 0.0, "p95": 0.0, "p99": 0.0, "max": 0.0, "mean": 0.0}
    return {
        "n": float(lengths.size),
        "min": float(lengths.min()),
        "p50": float(np.percentile(lengths, 50)),
        "p90": float(np.percentile(lengths, 90)),
        "p95": float(np.percentile(lengths, 95)),
        "p99": float(np.percentile(lengths, 99)),
        "max": float(lengths.max()),
        "mean": float(lengths.mean()),
    }


def make_folds(n: int, *, k: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    fold_id = np.empty(n, dtype=np.int16)
    for i, idx in enumerate(perm.tolist()):
        fold_id[idx] = i % k
    return fold_id


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=DEFAULT_FASTA, help="DM3Loc FASTA file")
    ap.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Output directory")
    ap.add_argument("--seed", type=int, default=42, help="Seed for deterministic CV folds")
    ap.add_argument("--folds", type=int, default=5, help="Number of CV folds (DM3Loc/LOCAS use 5)")
    ap.add_argument("--max-len", type=int, default=8000, help="Max sequence length after truncation")
    ap.add_argument("--left", type=int, default=4000, help="Left end length for truncation")
    ap.add_argument("--right", type=int, default=4000, help="Right end length for truncation")
    ap.add_argument(
        "--keep-cytoplasm-only",
        action="store_true",
        help="Keep sequences that only have the cytoplasm auxiliary label (no positives among 6 output classes)",
    )
    args = ap.parse_args()

    in_path: Path = args.input
    out_dir: Path = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if not in_path.exists():
        raise SystemExit(f"Input FASTA not found: {in_path}")

    records_path = out_dir / "dm3loc_trunc8000.jsonl.gz"
    folds_path = out_dir / "folds.json"
    summary_path = out_dir / "summary.json"
    readme_path = out_dir / "README.md"

    ids: List[str] = []
    accnums: List[str] = []
    lengths_orig: List[int] = []
    lengths_used: List[int] = []
    truncated_count = 0
    dropped_cytoplasm_only = 0
    labels6_list: List[List[int]] = []

    with gzip.open(records_path, "wt", encoding="utf-8") as out_f:
        for header, seq in iter_fasta(in_path):
            meta = parse_dm3loc_header(header)
            bits7 = meta["label_bits_7"]  # type: ignore[assignment]
            assert isinstance(bits7, list) and len(bits7) == 7

            labels6 = [int(bits7[i]) for i in IDX_6]
            cyt_aux = int(bits7[IDX_CYTOPLASM])

            if (not args.keep_cytoplasm_only) and sum(labels6) == 0:
                dropped_cytoplasm_only += 1
                continue

            seq_u, trunc = dm3loc_truncate(seq, max_len=args.max_len, left=args.left, right=args.right)
            if trunc != "none":
                truncated_count += 1

            acc = str(meta.get("ACCNUM", "")) if isinstance(meta.get("ACCNUM", ""), str) else ""
            gene_id = str(meta.get("Gene_ID", "")) if isinstance(meta.get("Gene_ID", ""), str) else ""
            gene_name = str(meta.get("Gene_Name", "")) if isinstance(meta.get("Gene_Name", ""), str) else ""
            rec_id = acc or header

            rec = {
                "id": rec_id,
                "accnum": acc or None,
                "gene_id": gene_id or None,
                "gene_name": gene_name or None,
                "sequence": seq_u,
                "length_original": int(len(str(seq).strip())),
                "length_used": int(len(seq_u)),
                "truncation": trunc,
                "label_names_6": LABEL_NAMES_6,
                "labels_6": labels6,
                "label_cytoplasm_aux": cyt_aux,
                "header": header,
            }
            out_f.write(json.dumps(rec) + "\n")

            ids.append(rec_id)
            accnums.append(acc or "")
            lengths_orig.append(rec["length_original"])
            lengths_used.append(rec["length_used"])
            labels6_list.append(labels6)

    n = len(ids)
    if n == 0:
        raise SystemExit("No records written; check input format and filters")

    fold_id = make_folds(n, k=int(args.folds), seed=int(args.seed))
    folds = {
        "k_folds": int(args.folds),
        "seed": int(args.seed),
        "note": "Deterministic fold assignment for 5-fold CV (DM3Loc/LOCAS-style). Fold i is the held-out test fold.",
        "folds": [
            {"fold": int(i), "test_ids": [ids[j] for j in range(n) if int(fold_id[j]) == i]}
            for i in range(int(args.folds))
        ],
    }
    folds_path.write_text(json.dumps(folds, indent=2), encoding="utf-8")

    y = np.asarray(labels6_list, dtype=np.int8)
    per_class = {LABEL_NAMES_6[i]: int(y[:, i].sum()) for i in range(y.shape[1])}
    summary = {
        "source": {
            "fasta": str(in_path),
            "format": "DM3Loc label bits in FASTA header",
            "label_order_7": ["Nucleus", "Exosome", "Cytosol", "Cytoplasm(aux)", "Ribosome", "Membrane", "ER"],
            "predicted_labels_6": LABEL_NAMES_6,
        },
        "preprocessing": {
            "max_len": int(args.max_len),
            "truncate": {"method": "end_concat", "left": int(args.left), "right": int(args.right)},
            "drop_cytoplasm_only": (not bool(args.keep_cytoplasm_only)),
        },
        "counts": {
            "n_sequences": int(n),
            "n_dropped_cytoplasm_only": int(dropped_cytoplasm_only),
            "n_truncated": int(truncated_count),
            "fraction_truncated": float(truncated_count) / float(n),
        },
        "lengths": {
            "original": compute_length_summary(np.asarray(lengths_orig, dtype=int)),
            "used": compute_length_summary(np.asarray(lengths_used, dtype=int)),
        },
        "labels": {
            "per_class_positive_counts": per_class,
            "num_labels_per_sample": compute_length_summary(y.sum(axis=1).astype(int)),
        },
        "outputs": {
            "records": str(records_path),
            "folds": str(folds_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    readme_path.write_text(
        "\n".join(
            [
                "# Task 2 — Subcellular Localization (DM3Loc/LOCAS benchmark)",
                "",
                "## Source dataset",
                f"- FASTA: `{in_path}`",
                "- Labels: 7-bit one-hot prefix in FASTA headers (DM3Loc convention).",
                "",
                "## Predicted labels (6)",
                "- " + ", ".join(LABEL_NAMES_6),
                "",
                "## Truncation (DM3Loc-style)",
                f"- If length > {args.max_len}: keep first {args.left} nt + last {args.right} nt (concatenate).",
                f"- Else: keep full sequence (no padding stored; pad in dataloader as needed).",
                "",
                "## Outputs",
                f"- Records: `{records_path}` (gzipped JSONL; one record per mRNA)",
                f"- Folds:   `{folds_path}` (deterministic {args.folds}-fold CV; seed={args.seed})",
                f"- Summary: `{summary_path}`",
                "",
                "## Notes",
                "- `label_cytoplasm_aux` is retained for reference only; cytoplasm is not part of the 6-label prediction target.",
                "- Some DM3Loc FASTA entries only have the cytoplasm auxiliary label; those are dropped by default (use `--keep-cytoplasm-only` to retain them).",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {n} records → {records_path}")
    print(f"Wrote folds       → {folds_path}")
    print(f"Wrote summary     → {summary_path}")
    print(f"Wrote README      → {readme_path}")


if __name__ == "__main__":
    main()
