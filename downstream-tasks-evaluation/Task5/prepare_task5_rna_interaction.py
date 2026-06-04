"""
Task 5: RNA-RNA Interaction Prediction — Dataset Preparation
=============================================================
Source:    DeepMirTar (Wen et al., Bioinformatics 2018)
           Used via RNAErnie's MirTarRAW benchmark dataset
Split:     RNAErnie convention (rna_rna_interaction.py):
             seed=0, 80% train, then 50/50 split of remainder → 80/10/10
Output:    Task5/splits/train.csv, val.csv, test.csv
           Task5/task5_summary.tsv
           Task5/task5_README.md
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
MIRTARRAW_CSV = Path(
    "/Users/andreyvasilyev/Documents/RNA-Protein-Research-Project-main"
    "/downstream-tasks-evaluation/Task5/data/MirTarRAW/MirTarRAW.csv"
)
OUTPUT_DIR = Path(
    "/Users/andreyvasilyev/Documents/RNA-Protein-Research-Project-main"
    "/downstream-tasks-evaluation/Task5"
)

# ---------------------------------------------------------------------------
# RNAErnie constants — from rna_rna_interaction.py source
# ---------------------------------------------------------------------------
SPLIT_RATIO = 0.80
SEED        = 0

# Padding lengths from miTAR paper (same dataset origin)
MIRNA_PAD_LEN = 26
MRNA_PAD_LEN  = 53
PAD_CHAR      = "N"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def pad_or_trim(seq: str, length: int) -> str:
    """Pad sequence with N to fixed length, or trim if longer."""
    seq = seq.strip().upper().replace("T", "U")
    if len(seq) >= length:
        return seq[:length]
    return seq + PAD_CHAR * (length - len(seq))


def main():
    if not MIRTARRAW_CSV.exists():
        print(f"ERROR: MirTarRAW.csv not found at:\n  {MIRTARRAW_CSV}")
        sys.exit(1)

    splits_dir = OUTPUT_DIR / "splits"
    splits_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Step 1: Load
    # ------------------------------------------------------------------ #
    print(f"Loading {MIRTARRAW_CSV.name} ...")
    df = pd.read_csv(MIRTARRAW_CSV)
    df.columns = [c.strip().lower() for c in df.columns]

    required = {"a_name", "a_seq", "b_name", "b_seq", "label"}
    missing_cols = required - set(df.columns)
    if missing_cols:
        print(f"ERROR: CSV is missing columns: {missing_cols}")
        print(f"Found columns: {list(df.columns)}")
        sys.exit(1)

    print(f"  Loaded {len(df):,} rows")
    print(f"  Positive pairs: {(df['label']==1).sum():,}")
    print(f"  Negative pairs: {(df['label']==0).sum():,}")

    # ------------------------------------------------------------------ #
    # Step 2: Pad sequences to fixed lengths (same as RNAErnie/miTAR)
    #         miRNA → 26 nt, mRNA site → 53 nt, T→U, pad with N
    # ------------------------------------------------------------------ #
    print(f"\n  miRNA length range before padding: "
          f"{df['a_seq'].str.len().min()}–{df['a_seq'].str.len().max()} nt")
    print(f"  mRNA  length range before padding: "
          f"{df['b_seq'].str.len().min()}–{df['b_seq'].str.len().max()} nt")

    df["a_seq"] = df["a_seq"].apply(lambda s: pad_or_trim(s, MIRNA_PAD_LEN))
    df["b_seq"] = df["b_seq"].apply(lambda s: pad_or_trim(s, MRNA_PAD_LEN))

    print(f"  After padding: miRNA={MIRNA_PAD_LEN} nt, mRNA={MRNA_PAD_LEN} nt")
    print(f"  Total rows kept: {len(df):,}")

    # ------------------------------------------------------------------ #
    # Step 3: Split — exact RNAErnie logic
    #   np_rng = np.random.RandomState(seed=0)
    #   np_rng.shuffle(data)
    #   train = data[:int(len*0.8)]
    #   test  = data[int(len*0.8):]
    #   We split test 50/50 into val/test → 80/10/10 total
    # ------------------------------------------------------------------ #
    data = df.values.tolist()
    rng  = np.random.RandomState(seed=SEED)
    rng.shuffle(data)

    n_total = len(data)
    n_train = int(n_total * SPLIT_RATIO)
    n_half  = (n_total - n_train) // 2

    cols = list(df.columns)
    train_df = pd.DataFrame(data[:n_train],                    columns=cols)
    val_df   = pd.DataFrame(data[n_train : n_train + n_half],  columns=cols)
    test_df  = pd.DataFrame(data[n_train + n_half :],          columns=cols)

    # ------------------------------------------------------------------ #
    # Step 4: Save
    # ------------------------------------------------------------------ #
    train_df.to_csv(splits_dir / "train.csv", index=False)
    val_df.to_csv(  splits_dir / "val.csv",   index=False)
    test_df.to_csv( splits_dir / "test.csv",  index=False)

    # ------------------------------------------------------------------ #
    # Step 5: Summary
    # ------------------------------------------------------------------ #
    print()
    rows = []
    for name, d in [("train", train_df), ("val", val_df), ("test", test_df)]:
        rows.append({
            "split":    name,
            "total":    len(d),
            "positive": int((d["label"] == 1).sum()),
            "negative": int((d["label"] == 0).sum()),
        })
        print(f"  {name:5s}  total={len(d):6,}  "
              f"pos={(d['label']==1).sum():6,}  "
              f"neg={(d['label']==0).sum():6,}")

    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "task5_summary.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ #
    # Step 6: README
    # ------------------------------------------------------------------ #
    with open(OUTPUT_DIR / "task5_README.md", "w") as f:
        f.write(f"""\
# Task 5 — RNA-RNA Interaction Prediction

## Source
DeepMirTar (Wen et al., Bioinformatics 2018)
Dataset: MirTarRAW (DeepMirTar + miRAW combined)
As used by: RNAErnie (Wang et al., Nature Machine Intelligence 2024)
  github.com/CatIIIIIIII/RNAErnie

## Task
Binary classification: predict whether a miRNA and an mRNA target site interact.

## Dataset stats
- Total pairs:               {n_total:,}
- Positive (interacting):    {(df['label']==1).sum():,}
- Negative (non-interacting):{(df['label']==0).sum():,}

## Sequence format
- miRNA (a_seq): padded/trimmed to {MIRNA_PAD_LEN} nt with N, T→U converted
- mRNA  (b_seq): padded/trimmed to {MRNA_PAD_LEN} nt with N, T→U converted
- label: 0 = no interaction, 1 = interaction

## Split — RNAErnie convention (rna_rna_interaction.py)
- Shuffle: numpy.random.RandomState(seed=0)
- Train: first 80%  → {len(train_df):,} pairs  (identical to RNAErnie train)
- Val:   next  10%  → {len(val_df):,} pairs
- Test:  last  10%  → {len(test_df):,} pairs   (comparable to RNAErnie results)

## Files
  splits/train.csv
  splits/val.csv
  splits/test.csv

## CSV columns
  a_name  miRNA identifier
  a_seq   miRNA sequence ({MIRNA_PAD_LEN} nt, N-padded, RNA)
  b_name  mRNA identifier
  b_seq   mRNA target site ({MRNA_PAD_LEN} nt, N-padded, RNA)
  label   0 or 1

## Load example
import pandas as pd
train = pd.read_csv("splits/train.csv")
X_mirna = train["a_seq"].values
X_mrna  = train["b_seq"].values
y       = train["label"].values
""")

    print(f"\n{'='*55}")
    print(f"Done!")
    print(f"  Splits  → {splits_dir}")
    print(f"  Summary → {OUTPUT_DIR / 'task5_summary.tsv'}")
    print(f"  README  → {OUTPUT_DIR / 'task5_README.md'}")
    print(f"{'='*55}")
    print()
    print("NOTE: Train is identical to RNAErnie's train set (seed=0, 80%).")
    print("      Test results are directly comparable to the RNAErnie paper.")


if __name__ == "__main__":
    main()