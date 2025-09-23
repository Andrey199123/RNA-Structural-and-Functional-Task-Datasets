#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import torch


def load_ids_from_npz(npz_path):
    arr = np.load(npz_path)
    ids = list(arr.files)
    return ids


def build_edge_index(train_csv, rna_ids, prot_ids):
    df = pd.read_csv(train_csv)
    rmap = {k: i for i, k in enumerate(rna_ids)}
    pmap = {k: i for i, k in enumerate(prot_ids)}
    src, dst = [], []
    offset = len(rna_ids)
    for rna, prot in zip(df['rna_id'], df['prot_id']):
        if rna in rmap and prot in pmap:
            u = rmap[rna]
            v = offset + pmap[prot]
            src.append(u); dst.append(v)
            src.append(v); dst.append(u)
    edge_index = torch.tensor([src, dst], dtype=torch.long)
    return edge_index


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_csv", required=True)
    ap.add_argument("--rna_npz", required=True)
    ap.add_argument("--prot_npz", required=True)
    ap.add_argument("--out_pt", required=True)
    args = ap.parse_args()

    rna_ids = load_ids_from_npz(args.rna_npz)
    prot_ids = load_ids_from_npz(args.prot_npz)
    edge_index = build_edge_index(args.train_csv, rna_ids, prot_ids)

    out_dir = Path(args.out_pt).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    torch.save({"edge_index": edge_index,
                "rna_ids": rna_ids,
                "prot_ids": prot_ids}, args.out_pt)
    print(f"Saved graph to {args.out_pt}")


if __name__ == "__main__":
    main()
