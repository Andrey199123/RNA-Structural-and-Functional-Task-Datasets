import sys
from pathlib import Path
import numpy as np
import pandas as pd
import h5py

THIS_DIR = Path(__file__).resolve().parent
CLIP_DIR = THIS_DIR / "data/PrismNet/data/clip_data"
OUTPUT_DIR = THIS_DIR

# 17 HeLa RBPs
HELA_RBPS = [
    "HNRNPC",
    "HNRNPU",
    "CSTF2",
    "METTL14",
    "METTL3",
    "TIA1",
    "TIAL1",
    "U2AF2",
    "U2AF65",
    "UPF1",
    "WTAP",
    "YTHDF2",
    "EIF4A3",
    "ELAVL1",
    "PTBP1",
    "PTBP1PTBP2",
    "eIF4AIII",
]

TRAIN_RATIO = 0.80
RANDOM_SEED = 42
SEQ_LEN     = 101
NUC_TO_IDX  = {"A": 0, "C": 1, "G": 2, "U": 3, "T": 3, "N": 4}

def parse_tsv(tsv_path: Path) -> pd.DataFrame:
    """
    PrismNet clip_data TSV columns:
      Type, Name, Seq, icshape, Score, label
    """
    df = pd.read_csv(tsv_path, sep="\t", usecols=["Seq", "icshape", "label"])
    df = df.rename(columns={"Seq": "sequence", "icshape": "structure"})
    df["sequence"] = df["sequence"].astype(str).str.upper().str.replace("T", "U", regex=False)
    df["structure"] = df["structure"].astype(str)
    df["label"] = pd.to_numeric(df["label"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["label"]).copy()
    df["label"] = df["label"].astype(np.int8)
    return df[["sequence", "structure", "label"]]


# Split — ERNIE-RNA: 80 / 10 / 10
def split_df(df: pd.DataFrame):
    rng     = np.random.default_rng(RANDOM_SEED)
    idx     = rng.permutation(len(df))
    n_train = int(round(TRAIN_RATIO * len(df)))
    n_half  = (len(df) - n_train) // 2
    return (
        df.iloc[idx[:n_train]].reset_index(drop=True),
        df.iloc[idx[n_train : n_train + n_half]].reset_index(drop=True),
        df.iloc[idx[n_train + n_half :]].reset_index(drop=True),
    )

# Encoding
def one_hot(seq: str) -> np.ndarray:
    arr = np.zeros((len(seq), 5), dtype=np.float32)
    for i, ch in enumerate(seq):
        arr[i, NUC_TO_IDX.get(ch, 4)] = 1.0
    return arr

def icshape_encode(struct: str) -> np.ndarray:
    # PrismNet stores icSHAPE as comma-separated floats (length 101); -1 means unmeasured.
    values = str(struct).split(",")
    arr = np.full(SEQ_LEN, -1.0, dtype=np.float32)
    for i, v in enumerate(values[:SEQ_LEN]):
        try:
            arr[i] = float(v)
        except ValueError:
            arr[i] = -1.0
    return arr

# Write HDF5
def write_h5(path: Path, train, val, test, rbp_name: str):
    with h5py.File(path, "w") as hf:
        hf.attrs["rbp_name"]     = rbp_name
        hf.attrs["cell_line"]    = "HeLa"
        hf.attrs["seq_length"]   = SEQ_LEN
        hf.attrs["split_method"] = "ERNIE-RNA 80/10/10 seed=42"
        hf.attrs["source"]       = "PrismNet (Xu et al., NAR 2023)"

        for name, df in [("train", train), ("val", val), ("test", test)]:
            g      = hf.create_group(name)
            seqs   = np.array([s.encode("ascii") for s in df["sequence"]], dtype="S101")
            oh     = np.stack([one_hot(s)        for s in df["sequence"]])
            ics    = np.stack([icshape_encode(s) for s in df["structure"]])
            labels = df["label"].to_numpy(dtype=np.int8)

            g.create_dataset("sequences", data=seqs,   compression="gzip")
            g.create_dataset("one_hot",   data=oh,     compression="gzip")
            g.create_dataset("icshape",   data=ics,    compression="gzip")
            g.create_dataset("labels",    data=labels, compression="gzip")
            g.attrs["n_total"]    = len(df)
            g.attrs["n_positive"] = int(labels.sum())
            g.attrs["n_negative"] = int((labels == 0).sum())


def main():
    if not CLIP_DIR.exists():
        print(f"ERROR: clip_data folder not found at:\n  {CLIP_DIR}")
        sys.exit(1)

    h5_dir = OUTPUT_DIR / "hdf5"
    h5_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading from: {CLIP_DIR}")
    print(f"Writing to:   {OUTPUT_DIR}\n")

    rows = []
    for rbp in HELA_RBPS:
        tsv_path = CLIP_DIR / f"{rbp}_Hela.tsv"
        if not tsv_path.exists():
            print(f"  ✗ {rbp:20s} — file not found, skipping")
            continue

        df = parse_tsv(tsv_path)
        df = df[df["sequence"].str.len() == SEQ_LEN].copy()
        if len(df) == 0:
            print(f"  ✗ {rbp:20s} — no valid rows, skipping")
            continue

        train, val, test = split_df(df)

        h5_path = h5_dir / f"{rbp}_Hela.h5"
        write_h5(h5_path, train, val, test, rbp)

        print(
            f"  ✓ {rbp:20s}  total={len(df):6d}  "
            f"train={len(train):5d}  val={len(val):4d}  test={len(test):4d}"
        )
        rows.append({
            "rbp": rbp, "total": len(df),
            "train": len(train), "val": len(val), "test": len(test),
            "train_pos": int(train["label"].sum()), "train_neg": int((train["label"]==0).sum()),
            "val_pos":   int(val["label"].sum()),   "val_neg":   int((val["label"]==0).sum()),
            "test_pos":  int(test["label"].sum()),  "test_neg":  int((test["label"]==0).sum()),
        })

    # Summary TSV
    summary_path = OUTPUT_DIR / "task4_summary.tsv"
    pd.DataFrame(rows).to_csv(summary_path, sep="\t", index=False)

    # README
    readme_path = OUTPUT_DIR / "task4_README.md"
    with open(readme_path, "w") as f:
        f.write(f"""\
# Task 4 — RNA-Protein Binding Prediction

## Source
PrismNet (Xu et al., Nucleic Acids Research 2023)
https://github.com/kuixu/PrismNet

## Split (ERNIE-RNA convention)
- Train: 80%
- Val:   10%
- Test:  10%
- Random seed: {RANDOM_SEED}

## Cell line
HeLa only

## HDF5 file structure
Each file hdf5/<RBP>_Hela.h5 has three groups: train / val / test
Each group contains:
  sequences  (N,)         raw 101-nt sequences as bytes
  one_hot    (N, 101, 5)  one-hot encoded ACGUN
  icshape    (N, 101)     icSHAPE scores, -1 = unmeasured
  labels     (N,)         0 = non-binding, 1 = binding

## Load example
import h5py
with h5py.File("hdf5/TIA1_Hela.h5", "r") as f:
    X = f["test/one_hot"][:]   # shape (N, 101, 5)
    y = f["test/labels"][:]    # shape (N,)

## RBPs processed ({len(rows)})
{", ".join(r["rbp"] for r in rows)}
""")

    print(f"\n{'='*55}")
    print(f"Done! {len(rows)}/17 RBPs processed.")
    print(f"HDF5 files → {h5_dir}")
    print(f"Summary    → {summary_path}")
    print(f"README     → {readme_path}")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
