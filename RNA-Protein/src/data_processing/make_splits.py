# This is from Phase 1 (ZHMolGraph Replication, not part of RNA MSA pipeline)
import argparse
import random
from pathlib import Path
import pandas as pd


def _detect_columns(df: pd.DataFrame):
    cols = {c.lower().strip(): c for c in df.columns}
    rna_keys  = ["rna_id", "rna", "rna names", "rna name"]
    prot_keys = ["prot_id", "protein", "protein names", "protein name"]
    rna_col = next((cols[k] for k in rna_keys if k in cols), None)
    prot_col = next((cols[k] for k in prot_keys if k in cols), None)
    if rna_col is None or prot_col is None:
        raise ValueError(f"Could not find RNA/protein columns in CSV. Got columns: {list(df.columns)}")
    return rna_col, prot_col


def random_edge_split(df, train_frac=0.8, val_frac=0.1, seed=42):
    df = df.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    n = len(df)
    n_train = int(train_frac * n)
    n_val = int(val_frac * n)
    return df.iloc[:n_train], df.iloc[n_train:n_train + n_val], df.iloc[n_train + n_val:]


def unknown_entity_split(df, entity_col, train_frac=0.8, val_frac=0.1, seed=42):
    random.seed(seed)
    entities = df[entity_col].unique().tolist()
    random.shuffle(entities)
    n = len(entities)
    n_train = int(train_frac * n)
    n_val = int(val_frac * n)
    train_e = set(entities[:n_train])
    val_e = set(entities[n_train:n_train + n_val])
    test_e = set(entities[n_train + n_val:])
    tr = df[df[entity_col].isin(train_e)]
    va = df[df[entity_col].isin(val_e)]
    te = df[df[entity_col].isin(test_e)]
    return tr, va, te


def unknown_nodes_split(df, seed=42):
    random.seed(seed)
    rnas = df['rna_id'].unique().tolist()
    prots = df['prot_id'].unique().tolist()
    random.shuffle(rnas)
    random.shuffle(prots)

    def thirds(xs):
        n = len(xs); a = n // 3; b = 2 * n // 3
        return set(xs[:a]), set(xs[a:b]), set(xs[b:])

    tr_r, va_r, te_r = thirds(rnas)
    tr_p, va_p, te_p = thirds(prots)
    tr = df[df['rna_id'].isin(tr_r) & df['prot_id'].isin(tr_p)]
    va = df[df['rna_id'].isin(va_r) & df['prot_id'].isin(va_p)]
    te = df[df['rna_id'].isin(te_r) & df['prot_id'].isin(te_p)]
    return tr, va, te


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs_csv", required=True, help="CSV with RNA/protein pairs (headers auto-detected)")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pairs_csv)
    rna_col, prot_col = _detect_columns(df)

    df = df[[rna_col, prot_col]].rename(columns={rna_col: "rna_id", prot_col: "prot_id"}).drop_duplicates()

    tr, va, te = random_edge_split(df, seed=args.seed)
    tr.to_csv(out_dir / "random_train.csv", index=False)
    va.to_csv(out_dir / "random_val.csv", index=False)
    te.to_csv(out_dir / "random_test.csv", index=False)

    tr, va, te = unknown_entity_split(df, "rna_id", seed=args.seed)
    tr.to_csv(out_dir / "uRNA_train.csv", index=False)
    va.to_csv(out_dir / "uRNA_val.csv", index=False)
    te.to_csv(out_dir / "uRNA_test.csv", index=False)

    tr, va, te = unknown_entity_split(df, "prot_id", seed=args.seed)
    tr.to_csv(out_dir / "uPROT_train.csv", index=False)
    va.to_csv(out_dir / "uPROT_val.csv", index=False)
    te.to_csv(out_dir / "uPROT_test.csv", index=False)

    tr, va, te = unknown_nodes_split(df, seed=args.seed)
    tr.to_csv(out_dir / "uNODES_train.csv", index=False)
    va.to_csv(out_dir / "uNODES_val.csv", index=False)
    te.to_csv(out_dir / "uNODES_test.csv", index=False)

    print(f"✅ Splits written to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
