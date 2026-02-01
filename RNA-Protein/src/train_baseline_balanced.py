#!/usr/bin/env python3
"""
Train ZHMolGraph baseline on the balanced dataset (5,916 pairs) for fair comparison
with the structural model.

This ensures apples-to-apples comparison since both models train on the same subset.
"""

import os
import sys
import math
import json
import random
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import torch
from torch import nn, optim
from tqdm import trange

THIS = Path(__file__).resolve()
SRC = THIS.parent
if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

from models.gnn_sage import RPISAGE
from models.vecnn_head import VecNNHead
from utils.io import load_npz_embeddings, set_seed

from sklearn.metrics import (
    matthews_corrcoef, accuracy_score, roc_auc_score, average_precision_score
)


def unknown_entity_split(df, entity_col, train_frac=0.8, val_frac=0.1, seed=42):
    """Split by holding out entities (RNAs or proteins)."""
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
    """Split where both RNA and protein are unseen in test."""
    random.seed(seed)
    rnas = df['rna_id'].unique().tolist()
    prots = df['prot_id'].unique().tolist()
    random.shuffle(rnas)
    random.shuffle(prots)

    def thirds(xs):
        n = len(xs)
        a = n // 3
        b = 2 * n // 3
        return set(xs[:a]), set(xs[a:b]), set(xs[b:])

    tr_r, va_r, te_r = thirds(rnas)
    tr_p, va_p, te_p = thirds(prots)
    tr = df[df['rna_id'].isin(tr_r) & df['prot_id'].isin(tr_p)]
    va = df[df['rna_id'].isin(va_r) & df['prot_id'].isin(va_p)]
    te = df[df['rna_id'].isin(te_r) & df['prot_id'].isin(te_p)]
    return tr, va, te


def build_edge_index(train_df, rna_ids, prot_ids):
    """Build edge index from training pairs."""
    rmap = {k: i for i, k in enumerate(rna_ids)}
    pmap = {k: i for i, k in enumerate(prot_ids)}
    src, dst = [], []
    offset = len(rna_ids)
    for rna, prot in zip(train_df['rna_id'], train_df['prot_id']):
        if rna in rmap and prot in pmap:
            u = rmap[rna]
            v = offset + pmap[prot]
            src.append(u)
            dst.append(v)
            src.append(v)
            dst.append(u)
    edge_index = torch.tensor([src, dst], dtype=torch.long)
    return edge_index


def pairs_to_index_tensors(df, rna_ids, prot_ids, device='cpu'):
    """Convert pairs dataframe to index tensors."""
    rmap = {k: i for i, k in enumerate(rna_ids)}
    pmap = {k: i for i, k in enumerate(prot_ids)}
    offset = len(rna_ids)
    u_list, v_list = [], []
    for rna, prot in zip(df['rna_id'], df['prot_id']):
        if rna in rmap and prot in pmap:
            u_list.append(rmap[rna])
            v_list.append(offset + pmap[prot])
    return torch.tensor(u_list, dtype=torch.long, device=device), \
           torch.tensor(v_list, dtype=torch.long, device=device)


def sample_negatives(num_rna, num_prot, num_samples, existing_set, seed=42):
    """Sample negative pairs."""
    rng = random.Random(seed)
    seen = set()
    out_u, out_v = [], []
    offset = num_rna
    max_trials = max(10_000, num_samples * 50)
    trials = 0
    while len(out_u) < num_samples and trials < max_trials:
        trials += 1
        u = rng.randrange(0, num_rna)
        v = offset + rng.randrange(0, num_prot)
        if (u, v) in existing_set or (u, v) in seen:
            continue
        seen.add((u, v))
        out_u.append(u)
        out_v.append(v)
    return torch.tensor(out_u, dtype=torch.long), torch.tensor(out_v, dtype=torch.long)


def build_existing_edge_set(df, rna_ids, prot_ids):
    """Build set of existing edges."""
    rmap = {k: i for i, k in enumerate(rna_ids)}
    pmap = {k: i for i, k in enumerate(prot_ids)}
    offset = len(rna_ids)
    s = set()
    for rna, prot in zip(df['rna_id'], df['prot_id']):
        if rna in rmap and prot in pmap:
            s.add((rmap[rna], offset + pmap[prot]))
    return s


def make_batches(pos_u, pos_v, batch_size):
    """Create batches from positive pairs."""
    n = pos_u.shape[0]
    if n == 0:
        return
    idx = torch.randperm(n)
    for i in range(0, n, batch_size):
        j = min(i + batch_size, n)
        yield pos_u[idx[i:j]], pos_v[idx[i:j]]


def _metrics_from_arrays(y_true, y_prob, thr=0.5):
    """Compute metrics from arrays."""
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    y_pred = (y_prob >= thr).astype(int)
    if y_true.size == 0:
        return {"acc": 0.0, "mcc": 0.0, "auroc": 0.0, "auprc": 0.0}
    return {
        "acc": accuracy_score(y_true, y_pred),
        "mcc": matthews_corrcoef(y_true, y_pred),
        "auroc": roc_auc_score(y_true, y_prob),
        "auprc": average_precision_score(y_true, y_prob),
    }


def _best_mcc_threshold(y_true, y_prob, grid=None):
    """Find threshold that maximizes MCC."""
    if grid is None:
        grid = np.linspace(0.05, 0.95, 37)
    best_t, best_m = 0.5, -1.0
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    for t in grid:
        pred = (y_prob >= t).astype(int)
        m = matthews_corrcoef(y_true, pred)
        if m > best_m:
            best_m, best_t = m, t
    return best_t, best_m


def evaluate(model, head, x_rna, x_prot, edge_index, df, rna_ids, prot_ids, device):
    """Evaluate model on a dataframe of pairs (using actual labels, not sampled negatives)."""
    model.eval()
    head.eval()
    
    with torch.no_grad():
        z = model(x_rna, x_prot, edge_index)
        
        # Get indices for all pairs
        rmap = {k: i for i, k in enumerate(rna_ids)}
        pmap = {k: i for i, k in enumerate(prot_ids)}
        offset = len(rna_ids)
        
        u_list, v_list, labels = [], [], []
        for _, row in df.iterrows():
            rna = row['rna_id']
            prot = row['prot_id']
            label = row['Labels']
            if rna in rmap and prot in pmap:
                u_list.append(rmap[rna])
                v_list.append(offset + pmap[prot])
                labels.append(label)
        
        u = torch.tensor(u_list, dtype=torch.long, device=device)
        v = torch.tensor(v_list, dtype=torch.long, device=device)
        
        scores = torch.sigmoid(head(z[u], z[v])).cpu().numpy()
    
    y_true = np.array(labels)
    y_prob = scores
    
    return {"y_true": y_true, "y_prob": y_prob}


def main():
    # Get script directory for relative paths
    script_dir = Path(__file__).parent.parent  # RNA-Protein/
    
    parser = argparse.ArgumentParser(description="Train ZHMolGraph baseline on balanced dataset")
    parser.add_argument('--data-dir', type=Path, default=None)
    parser.add_argument('--split', type=str, default='random',
                       choices=['random', 'uRNA', 'uPROT', 'uNODES'])
    parser.add_argument('--epochs', type=int, default=30)
    parser.add_argument('--batch-size', type=int, default=1024)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--hidden-dim', type=int, default=128)
    parser.add_argument('--out-dim', type=int, default=64)
    parser.add_argument('--num-layers', type=int, default=2)
    parser.add_argument('--dropout', type=float, default=0.2)
    parser.add_argument('--head-hidden', type=int, default=256)
    parser.add_argument('--neg-ratio', type=float, default=1.0)
    parser.add_argument('--patience', type=int, default=5)
    parser.add_argument('--device', type=str, default='auto')
    parser.add_argument('--save-dir', type=Path, default=None)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    
    # Set default paths
    if args.data_dir is None:
        args.data_dir = script_dir / 'data'
    if args.save_dir is None:
        args.save_dir = script_dir / 'checkpoints_baseline_balanced'
    
    set_seed(args.seed)
    
    # Device
    if args.device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(args.device)
    print(f"Using device: {device}")
    
    # Paths
    data_dir = args.data_dir
    pairs_csv = data_dir / 'RPI7317/rpi7317_pairs_balanced.csv'
    rna_npz = data_dir / 'embeddings/rna_fm_rpi7317.npz'
    prot_npz = data_dir / 'embeddings/esm2_35m_rpi7317.npz'
    
    # Check files exist
    if not pairs_csv.exists():
        print(f"ERROR: Balanced pairs file not found: {pairs_csv}")
        print("Please run generate_negative_samples.py first.")
        return
    
    # Load embeddings
    print("Loading embeddings...")
    (rna_ids, Xr), (prot_ids, Xp) = load_npz_embeddings(str(rna_npz), str(prot_npz))
    x_rna = Xr.to(device)
    x_prot = Xp.to(device)
    print(f"  RNA embeddings: {len(rna_ids)} x {x_rna.shape[1]}")
    print(f"  Protein embeddings: {len(prot_ids)} x {x_prot.shape[1]}")
    
    # Load balanced pairs
    print("Loading balanced pairs...")
    full_df = pd.read_csv(pairs_csv, encoding='utf-8-sig')
    full_df = full_df.rename(columns={'RNA names': 'rna_id', 'Protein names': 'prot_id'})
    print(f"  Total pairs: {len(full_df)}")
    print(f"  Positive: {(full_df['Labels'] == 1).sum()}")
    print(f"  Negative: {(full_df['Labels'] == 0).sum()}")

    
    # Create splits
    print(f"Creating {args.split} split...")
    if args.split == 'random':
        full_df = full_df.sample(frac=1, random_state=args.seed).reset_index(drop=True)
        n = len(full_df)
        train_df = full_df[:int(0.8*n)]
        val_df = full_df[int(0.8*n):int(0.9*n)]
        test_df = full_df[int(0.9*n):]
    elif args.split == 'uRNA':
        train_df, val_df, test_df = unknown_entity_split(full_df, 'rna_id', seed=args.seed)
    elif args.split == 'uPROT':
        train_df, val_df, test_df = unknown_entity_split(full_df, 'prot_id', seed=args.seed)
    elif args.split == 'uNODES':
        train_df, val_df, test_df = unknown_nodes_split(full_df, seed=args.seed)
    
    print(f"  Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")
    
    # Filter to only positive pairs for GNN training (negatives sampled during training)
    train_pos = train_df[train_df['Labels'] == 1].copy()
    print(f"  Training positive pairs: {len(train_pos)}")
    
    # Build edge index from positive training pairs
    edge_index = build_edge_index(train_pos, rna_ids, prot_ids).to(device)
    print(f"  Edge index shape: {edge_index.shape}")
    
    # Get positive pair indices for training
    pos_u_tr, pos_v_tr = pairs_to_index_tensors(train_pos, rna_ids, prot_ids, device=device)
    
    # Build existing edge set (all positive pairs)
    all_pos = full_df[full_df['Labels'] == 1]
    exist_all = build_existing_edge_set(all_pos, rna_ids, prot_ids)
    
    # Create model
    print("Creating model...")
    model = RPISAGE(
        in_rna=x_rna.shape[1],
        in_prot=x_prot.shape[1],
        hidden_dim=args.hidden_dim,
        out_dim=args.out_dim,
        num_layers=args.num_layers,
        dropout=args.dropout,
    ).to(device)
    head = VecNNHead(args.out_dim, hidden=args.head_hidden).to(device)
    
    total_params = sum(p.numel() for p in model.parameters()) + sum(p.numel() for p in head.parameters())
    print(f"  Model parameters: {total_params:,}")
    
    # Training setup
    params = list(model.parameters()) + list(head.parameters())
    opt = optim.Adam(params, lr=args.lr, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()
    
    args.save_dir.mkdir(parents=True, exist_ok=True)
    print(f"Checkpoints will be saved to: {args.save_dir.absolute()}")
    
    best_val_mcc = -1.0
    best_snapshot_path = args.save_dir / f'best_model_{args.split}.pt'
    bad_epochs = 0
    
    # Training loop
    print("\nTraining...")
    for epoch in trange(1, args.epochs + 1, desc="epochs"):
        model.train()
        head.train()
        total_loss = 0.0
        
        for pu, pv in make_batches(pos_u_tr, pos_v_tr, args.batch_size):
            num_neg = int(math.ceil(len(pu) * args.neg_ratio))
            neg_u, neg_v = sample_negatives(
                x_rna.shape[0], x_prot.shape[0], num_neg, exist_all,
                seed=random.randint(0, 10**9)
            )
            pu = pu.to(device)
            pv = pv.to(device)
            neg_u = neg_u.to(device)
            neg_v = neg_v.to(device)
            
            z = model(x_rna, x_prot, edge_index)
            pos_logit = head(z[pu], z[pv])
            neg_logit = head(z[neg_u], z[neg_v])
            
            logits = torch.cat([pos_logit, neg_logit], dim=0)
            labels = torch.cat([torch.ones_like(pos_logit), torch.zeros_like(neg_logit)], dim=0)
            
            loss = loss_fn(logits, labels)
            opt.zero_grad()
            loss.backward()
            opt.step()
            total_loss += loss.item()
        
        # Evaluate on validation and test sets (using actual labels)
        val_raw = evaluate(model, head, x_rna, x_prot, edge_index, val_df, rna_ids, prot_ids, device)
        test_raw = evaluate(model, head, x_rna, x_prot, edge_index, test_df, rna_ids, prot_ids, device)
        
        thr, _ = _best_mcc_threshold(val_raw["y_true"], val_raw["y_prob"])
        val_m = _metrics_from_arrays(val_raw["y_true"], val_raw["y_prob"], thr=thr)
        test_m = _metrics_from_arrays(test_raw["y_true"], test_raw["y_prob"], thr=thr)
        
        print(f"\nEpoch {epoch} | loss={total_loss:.4f} "
              f"| VAL acc={val_m['acc']:.4f} mcc={val_m['mcc']:.4f} auroc={val_m['auroc']:.4f} "
              f"| TEST acc={test_m['acc']:.4f} mcc={test_m['mcc']:.4f} auroc={test_m['auroc']:.4f} "
              f"| thr={thr:.2f}")
        
        if val_m["mcc"] > best_val_mcc:
            best_val_mcc = val_m["mcc"]
            bad_epochs = 0
            torch.save({
                "model": model.state_dict(),
                "head": head.state_dict(),
                "args": vars(args),
                "rna_ids": rna_ids,
                "prot_ids": prot_ids,
                "best_thr": float(thr),
                "best_val_metrics": val_m,
                "test_metrics": test_m,
            }, best_snapshot_path)
        else:
            bad_epochs += 1
            if bad_epochs >= args.patience:
                print(f"Early stopping (no val MCC improvement for {args.patience} epochs).")
                break
    
    # Load best model and evaluate
    print("\nLoading best model for final evaluation...")
    checkpoint = torch.load(best_snapshot_path, weights_only=False)
    model.load_state_dict(checkpoint["model"])
    head.load_state_dict(checkpoint["head"])
    best_thr = checkpoint["best_thr"]
    
    test_raw = evaluate(model, head, x_rna, x_prot, edge_index, test_df, rna_ids, prot_ids, device)
    test_m = _metrics_from_arrays(test_raw["y_true"], test_raw["y_prob"], thr=best_thr)
    
    print("\n" + "="*60)
    print("FINAL TEST RESULTS (ZHMolGraph Baseline on Balanced Dataset)")
    print("="*60)
    print(f"Split: {args.split}")
    print(f"Accuracy: {test_m['acc']:.4f}")
    print(f"MCC:      {test_m['mcc']:.4f}")
    print(f"AUROC:    {test_m['auroc']:.4f}")
    print(f"AUPRC:    {test_m['auprc']:.4f}")
    print(f"Threshold: {best_thr:.2f}")
    
    # Save results
    results = {
        'split': args.split,
        'test_metrics': test_m,
        'best_val_mcc': best_val_mcc,
        'best_threshold': best_thr,
        'dataset': 'rpi7317_balanced',
        'n_pairs': len(full_df),
        'n_train': len(train_df),
        'n_val': len(val_df),
        'n_test': len(test_df),
    }
    results_path = args.save_dir / f'results_{args.split}.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {results_path}")


if __name__ == "__main__":
    main()
