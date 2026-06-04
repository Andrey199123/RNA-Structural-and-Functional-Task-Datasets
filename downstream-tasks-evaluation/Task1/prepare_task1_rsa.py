import json
import tarfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / "preprocessed"

RNASOL_DIR = SCRIPT_DIR / "data" / "rnasol_benchmark"
RNASOL_EXTRACTED_DIR = RNASOL_DIR / "extracted"

TR89_TINKER_FILE = RNASOL_EXTRACTED_DIR / "datasets" / "tr89" / "tr89_tinker_data.txt"
TS44_TINKER_FILE = RNASOL_EXTRACTED_DIR / "ts44" / "ts44_tinker_data.txt"
CN48_TINKER_FILE = RNASOL_EXTRACTED_DIR / "cn48" / "cn48_tinker_data.txt"

BASE_MAX_ASA = {"A": 400.0, "G": 400.0, "U": 350.0, "C": 350.0}


def normalize_asa_with_mask(asa_values: List[float], sequence: str) -> Tuple[List[float], List[int]]:
    if len(asa_values) != len(sequence):
        raise ValueError(f"Length mismatch: seq={len(sequence)} asa={len(asa_values)}")

    norm: List[float] = []
    mask: List[int] = []
    for base, asa in zip(sequence, asa_values):
        if asa < 0:
            norm.append(0.0)
            mask.append(0)
            continue
        scale = BASE_MAX_ASA.get(base.upper(), 400.0)
        norm.append(min(float(asa) / float(scale), 1.0))
        mask.append(1)
    return norm, mask


def _parse_rsa_flags(s: str) -> List[int]:
    out: List[int] = []
    i = 0
    while i < len(s):
        if s[i] == "-" and i + 1 < len(s) and s[i + 1] == "1":
            out.append(-1)
            i += 2
        else:
            ch = s[i]
            if ch not in ("0", "1"):
                raise ValueError(f"Unexpected RSA flag character at {i}: {s[i:i+10]!r}")
            out.append(int(ch))
            i += 1
    return out


def parse_rnasol_benchmark_txt(text: str) -> List[Dict]:
    records: List[Dict] = []
    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    i = 0

    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue
        if not line.startswith(">"):
            raise ValueError(f"Unexpected line (expected '>id'): {line[:60]!r}")
        if i + 3 >= len(lines):
            raise ValueError(f"Truncated record at line {i}")

        rec_id = line[1:].strip()
        seq = lines[i + 1].strip().upper().replace("T", "U")
        asa_raw = [float(x) for x in lines[i + 2].strip().split()]
        rsa_flags = _parse_rsa_flags(lines[i + 3].strip())

        if len(seq) != len(asa_raw):
            raise ValueError(f"Length mismatch for {rec_id}: seq={len(seq)} asa={len(asa_raw)}")
        if len(seq) != len(rsa_flags):
            raise ValueError(f"Length mismatch for {rec_id}: seq={len(seq)} rsa_flags={len(rsa_flags)}")

        asa_norm, mask = normalize_asa_with_mask(asa_raw, seq)

        for j, (a, rf, m) in enumerate(zip(asa_raw, rsa_flags, mask)):
            if (a < 0) != (rf == -1) or (a < 0) != (m == 0):
                raise ValueError(
                    f"Inconsistent missing label for {rec_id} pos {j}: asa={a} rsa_flag={rf} mask={m}"
                )

        records.append(
            {
                "id": rec_id,
                "full_id": rec_id,
                "sequence": seq,
                "asa_raw": asa_raw,
                "asa_normalized": asa_norm,
                "label_mask": mask,
                "length": len(seq),
            }
        )
        i += 4

    return records


def _ensure_rnasol_extracted() -> None:
    if TR89_TINKER_FILE.exists() and TS44_TINKER_FILE.exists() and CN48_TINKER_FILE.exists():
        return

    tarballs = [
        RNASOL_DIR / "tr89.tar.bz2",
        RNASOL_DIR / "ts44.tar.bz2",
        RNASOL_DIR / "cn48.tar.bz2",
    ]
    missing = [str(t) for t in tarballs if not t.exists()]
    if missing:
        raise FileNotFoundError("Missing RNAsol tarballs:\n  " + "\n  ".join(missing))

    RNASOL_EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)
    for t in tarballs:
        with tarfile.open(t, "r:*") as tf:
            try:
                tf.extractall(path=RNASOL_EXTRACTED_DIR, filter="data")
            except TypeError:
                tf.extractall(path=RNASOL_EXTRACTED_DIR)


def process_dataset() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    _ensure_rnasol_extracted()

    splits = {
        "train": parse_rnasol_benchmark_txt(TR89_TINKER_FILE.read_text(encoding="utf-8", errors="replace")),
        "val": [],
        "test_ts44": parse_rnasol_benchmark_txt(TS44_TINKER_FILE.read_text(encoding="utf-8", errors="replace")),
        "test_cn48": parse_rnasol_benchmark_txt(CN48_TINKER_FILE.read_text(encoding="utf-8", errors="replace")),
    }

    dataset: Dict[str, List[Dict]] = {}
    statistics: Dict[str, Dict] = {}

    for split_name, records in splits.items():
        asa_stats = {"min": float("inf"), "max": float("-inf"), "mean": 0.0, "count": 0}
        seq_length_stats: List[int] = []
        missing_count = 0

        for rec in records:
            asa_values = rec["asa_raw"]
            mask = rec["label_mask"]
            for a, m in zip(asa_values, mask):
                if not m:
                    missing_count += 1
                    continue
                asa_stats["min"] = min(asa_stats["min"], a)
                asa_stats["max"] = max(asa_stats["max"], a)
                asa_stats["mean"] += a
                asa_stats["count"] += 1
            seq_length_stats.append(len(rec["sequence"]))

        if asa_stats["count"] > 0:
            asa_stats["mean"] /= asa_stats["count"]
        else:
            asa_stats = {"min": 0.0, "max": 0.0, "mean": 0.0, "count": 0}

        statistics[split_name] = {
            "num_sequences": len(records),
            "num_missing_labels": missing_count,
            "asa_stats": asa_stats,
            "seq_length_stats": {
                "min": int(min(seq_length_stats)) if seq_length_stats else 0,
                "max": int(max(seq_length_stats)) if seq_length_stats else 0,
                "mean": float(np.mean(seq_length_stats)) if seq_length_stats else 0.0,
                "median": float(np.median(seq_length_stats)) if seq_length_stats else 0.0,
            },
        }

        dataset[split_name] = records

    for split_name, split_data in dataset.items():
        output_file = OUTPUT_DIR / f"{split_name}.json"
        with open(output_file, "w") as f:
            json.dump(split_data, f, indent=2)

    stats_file = OUTPUT_DIR / "statistics.json"
    with open(stats_file, "w") as f:
        json.dump(statistics, f, indent=2)

if __name__ == "__main__":
    process_dataset()