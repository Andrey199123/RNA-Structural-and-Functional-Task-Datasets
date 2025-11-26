# This is from Phase 1 (ZHMolGraph Replication, not part of RNA MSA pipeline)
import os, sys, json, random
from pathlib import Path
import argparse
import yaml
import numpy as np
import pandas as pd
import torch
from tqdm import trange
from sklearn.model_selection import KFold

THIS = Path(__file__).resolve()
SRC = THIS.parents[1]
if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

from models.gnn_sage import RPISAGE
from models.vecnn_head import VecNNHead
from utils.io import load_npz_embeddings, pairs_to_index_tensors, set_seed
from train_linkpred import (
    build_edge_index, sample_negatives, build_existing_edge_set,
    evaluate, _best_mcc_threshold, _metrics_from_arrays
)

from sklearn.metrics import (
    matthews_corrcoef, accuracy_score, roc_auc_score, average_precision_score
)

def device_from_cfg(cfg):
    dev = cfg.get("runtime", {}).get("device", "auto")
    if dev == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(dev)

def run_single_fold(fold_idx, train_df, val_df, test_df, rna_ids, prot_ids, 
                   x_rna, x_prot, cfg, device):
    set_seed(cfg["runtime"]["seed"] + fold_idx)
    
    edge_index = build_edge_index(train_df, rna_ids, prot_ids).to(device)
    
    pos_u_tr, pos_v_tr = pairs_to_index_tensors(train_df, rna_ids, prot_ids, device=device)
    
    all_df = pd.concat([train_df, val_df, test_df], ignore_index=True)
    exist_all = build_existing_edge_set(all_df, rna_ids, prot_ids)
    
    model = RPISAGE(
        in_rna=x_rna.shape[1],
        in_prot=x_prot.shape[1],
        hidden_dim=cfg["model"]["hidden_dim"],
        out_dim=cfg["model"]["out_dim"],
        num_layers=cfg["model"]["num_layers"],
        dropout=cfg["model"]["dropout"],
    ).to(device)

    head_hidden = cfg["model"].get("head_hidden", 256)
    head = VecNNHead(cfg["model"]["out_dim"], hidden=head_hidden).to(device)
    
    params = list(model.parameters()) + list(head.parameters())
    opt = torch.optim.Adam(params, lr=cfg["train"]["lr"], weight_decay=cfg["train"]["weight_decay"])
    loss_fn = torch.nn.BCEWithLogitsLoss()
    
    best_val_mcc = -1.0
    patience = cfg["train"]["early_stop_patience"]
    bad_epochs = 0
    
    epochs = cfg["train"]["epochs"]
    bsz = cfg["train"]["batch_size"]
    neg_ratio = cfg["train"]["neg_ratio"]
    eval_neg_ratio = cfg["train"]["eval_neg_ratio"]
    
    for epoch in range(1, epochs + 1):
        model.train(); head.train()
        total_loss = 0.0
        
        n = pos_u_tr.shape[0]
        if n == 0:
            break
        idx = torch.randperm(n)

        for i in range(0, n, bsz):

            j = min(i + bsz, n)
            pu = pos_u_tr[idx[i:j]]
            pv = pos_v_tr[idx[i:j]]
            
            num_neg = int(torch.ceil(torch.tensor(len(pu) * neg_ratio)).item())
            neg_u, neg_v = sample_negatives(
                x_rna.shape[0], x_prot.shape[0], num_neg, exist_all,
                seed=random.randint(0, 10**9)
            )
            pu = pu.to(device); pv = pv.to(device)
            neg_u = neg_u.to(device); neg_v = neg_v.to(device)
            
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
        
        val_raw = evaluate(
            model, head, x_rna, x_prot, edge_index, val_df,
            rna_ids, prot_ids, eval_neg_ratio, exist_all, device
        )
        
        thr, best_mcc_val = _best_mcc_threshold(val_raw["y_true"], val_raw["y_prob"])
        val_m = _metrics_from_arrays(val_raw["y_true"], val_raw["y_prob"], thr=thr)
        
        if val_m["mcc"] > best_val_mcc:

            best_val_mcc = val_m["mcc"]
            bad_epochs = 0
            best_model_state = {k: v.cpu() for k, v in model.state_dict().items()}
            best_head_state = {k: v.cpu() for k, v in head.state_dict().items()}
            best_thr = thr

        else:
            bad_epochs += 1
            if bad_epochs >= patience:
                break
    
    model.load_state_dict({k: v.to(device) for k, v in best_model_state.items()})
    head.load_state_dict({k: v.to(device) for k, v in best_head_state.items()})
    
    test_raw = evaluate(
        model, head, x_rna, x_prot, edge_index, test_df,
        rna_ids, prot_ids, eval_neg_ratio, exist_all, device
    )
    test_m = _metrics_from_arrays(test_raw["y_true"], test_raw["y_prob"], thr=best_thr)
    
    return {
        "fold": fold_idx,
        "val_mcc": best_val_mcc,
        "test_acc": test_m["acc"],
        "test_mcc": test_m["mcc"],
        "test_auroc": test_m["auroc"],
        "test_auprc": test_m["auprc"],
        "threshold": best_thr
    }


def run_cross_validation(dataset_name, split_type, cfg, n_folds=5):

    print(f"\n{'='*60}")
    print(f"Running {n_folds}-fold CV on {dataset_name} with {split_type} split")
    print(f"{'='*60}")
    
    (rna_ids, Xr), (prot_ids, Xp) = load_npz_embeddings(
        cfg["embeddings"]["rna_npz"],
        cfg["embeddings"]["prot_npz"]
    )
    device = device_from_cfg(cfg)
    x_rna = Xr.to(device)
    x_prot = Xp.to(device)
    
    pairs_path = cfg["data"]["pairs_csv"]
    all_pairs = pd.read_csv(pairs_path)
    
    cols = {c.lower().strip(): c for c in all_pairs.columns}
    rna_keys = ["rna_id", "rna", "rna names", "rna name"]
    prot_keys = ["prot_id", "protein", "protein names", "protein name"]
    rna_col = next((cols[k] for k in rna_keys if k in cols), None)
    prot_col = next((cols[k] for k in prot_keys if k in cols), None)
    
    if rna_col is None or prot_col is None:
        raise ValueError(f"Could not find RNA/protein columns in {pairs_path}")
    
    df = all_pairs[[rna_col, prot_col]].rename(columns={rna_col: "rna_id", prot_col: "prot_id"}).drop_duplicates()
    
    rna_set = set(rna_ids)
    prot_set = set(prot_ids)
    df = df[df['rna_id'].isin(rna_set) & df['prot_id'].isin(prot_set)]
    
    print(f"Total pairs after filtering: {len(df)}")
    print(f"Unique RNAs: {df['rna_id'].nunique()}")
    print(f"Unique Proteins: {df['prot_id'].nunique()}")
    
    if split_type == "random":

        kf = KFold(n_splits=n_folds, shuffle=True, random_state=cfg["runtime"]["seed"])
        edge_splits = list(kf.split(df))
        splits = []
        for train_idx, test_idx in edge_splits:
            train_df = df.iloc[train_idx]
            test_df = df.iloc[test_idx]
            # Split train into train/val
            val_size = len(train_df) // 10
            val_df = train_df.sample(n=val_size, random_state=cfg["runtime"]["seed"])
            train_df = train_df.drop(val_df.index)
            splits.append((train_df.index, val_df.index, test_df.index))

    elif split_type == "uRNA":

        # Unknown RNA split
        rnas = df['rna_id'].unique()
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=cfg["runtime"]["seed"])
        rna_splits = list(kf.split(rnas))
        splits = []

        for train_rna_idx, test_rna_idx in rna_splits:

            train_rnas = set(rnas[train_rna_idx])
            test_rnas = set(rnas[test_rna_idx])
            train_df = df[df['rna_id'].isin(train_rnas)]
            test_df = df[df['rna_id'].isin(test_rnas)]

            # Split train into train/val
            val_size = len(train_df) // 10
            val_df = train_df.sample(n=val_size, random_state=cfg["runtime"]["seed"])
            train_df = train_df.drop(val_df.index)
            splits.append((train_df.index, val_df.index, test_df.index))

    elif split_type == "uPROT":

        # Unknown Protein split
        prots = df['prot_id'].unique()
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=cfg["runtime"]["seed"])
        prot_splits = list(kf.split(prots))
        splits = []

        for train_prot_idx, test_prot_idx in prot_splits:
            train_prots = set(prots[train_prot_idx])
            test_prots = set(prots[test_prot_idx])
            train_df = df[df['prot_id'].isin(train_prots)]
            test_df = df[df['prot_id'].isin(test_prots)]

            # Split train into train/val
            val_size = len(train_df) // 10
            val_df = train_df.sample(n=val_size, random_state=cfg["runtime"]["seed"])
            train_df = train_df.drop(val_df.index)
            splits.append((train_df.index, val_df.index, test_df.index))

    elif split_type == "uNODES":

        # Unknown Nodes split (both RNA and protein held out)
        rnas = df['rna_id'].unique()
        prots = df['prot_id'].unique()
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=cfg["runtime"]["seed"])
        rna_splits = list(kf.split(rnas))
        prot_splits = list(kf.split(prots))
        splits = []

        for (train_rna_idx, test_rna_idx), (train_prot_idx, test_prot_idx) in zip(rna_splits, prot_splits):
            
            train_rnas = set(rnas[train_rna_idx])
            test_rnas = set(rnas[test_rna_idx])
            train_prots = set(prots[train_prot_idx])
            test_prots = set(prots[test_prot_idx])
            train_df = df[df['rna_id'].isin(train_rnas) & df['prot_id'].isin(train_prots)]
            test_df = df[df['rna_id'].isin(test_rnas) & df['prot_id'].isin(test_prots)]
            # Split train into train/val
            val_size = len(train_df) // 10
            val_df = train_df.sample(n=val_size, random_state=cfg["runtime"]["seed"])
            train_df = train_df.drop(val_df.index)
            splits.append((train_df.index, val_df.index, test_df.index))

    else:
        raise ValueError(f"Unknown split type: {split_type}")
    
    # Run cross-validation
    results = []
    for fold_idx, (train_idx, val_idx, test_idx) in enumerate(splits):

        print(f"\nFold {fold_idx + 1}/{n_folds}")
        train_df = df.iloc[train_idx]
        val_df = df.iloc[val_idx]
        test_df = df.iloc[test_idx]
        
        print(f"Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")
        
        result = run_single_fold(fold_idx, train_df, val_df, test_df, 
                               rna_ids, prot_ids, x_rna, x_prot, cfg, device)
        results.append(result)
        print(f"Fold {fold_idx + 1} - Val MCC: {result['val_mcc']:.4f}, "
              f"Test MCC: {result['test_mcc']:.4f}, Test AUROC: {result['test_auroc']:.4f}")
    
    results_df = pd.DataFrame(results)
    summary = {
        "dataset": dataset_name,
        "split_type": split_type,
        "n_folds": n_folds,
        "mean_test_acc": results_df["test_acc"].mean(),
        "std_test_acc": results_df["test_acc"].std(),
        "mean_test_mcc": results_df["test_mcc"].mean(),
        "std_test_mcc": results_df["test_mcc"].std(),
        "mean_test_auroc": results_df["test_auroc"].mean(),
        "std_test_auroc": results_df["test_auroc"].std(),
        "mean_test_auprc": results_df["test_auprc"].mean(),
        "std_test_auprc": results_df["test_auprc"].std(),
        "fold_results": results
    }
    
    print(f"\n{'-'*40}")
    print(f"SUMMARY - {dataset_name} ({split_type})")
    print(f"{'-'*40}")
    print(f"Test Accuracy:  {summary['mean_test_acc']:.4f} ± {summary['std_test_acc']:.4f}")
    print(f"Test MCC:       {summary['mean_test_mcc']:.4f} ± {summary['std_test_mcc']:.4f}")
    print(f"Test AUROC:     {summary['mean_test_auroc']:.4f} ± {summary['std_test_auroc']:.4f}")
    print(f"Test AUPRC:     {summary['mean_test_auprc']:.4f} ± {summary['std_test_auprc']:.4f}")
    
    return summary


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["rpi7317", "npinter2"], required=True)
    ap.add_argument("--split_type", choices=["random", "uRNA", "uPROT", "uNODES"], required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", help="Output JSON file for results")
    args = ap.parse_args()
    
    cfg = yaml.safe_load(Path(args.config).read_text())
    
    # Set up paths based on dataset
    if args.dataset == "rpi7317":
        cfg["data"]["pairs_csv"] = "data/RPI7317/rpi7317_pairs.csv"
        cfg["embeddings"]["rna_npz"] = "data/embeddings/rna_fm_rpi7317.npz"
        cfg["embeddings"]["prot_npz"] = "data/embeddings/esm2_35m_rpi7317.npz"
    elif args.dataset == "npinter2":
        # Haven't used NPInter yet but it's here just in case
        cfg["data"]["pairs_csv"] = "data/raw/NPInter10412_dataset.txt"
        cfg["embeddings"]["rna_npz"] = "data/embeddings/rna_fm_npinter2.npz"
        cfg["embeddings"]["prot_npz"] = "data/embeddings/esm2_35m_npinter2.npz"
    
    # Run cross-validation
    summary = run_cross_validation(args.dataset, args.split_type, cfg)
    
    # Save results
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"\nResults saved to: {args.output}")
    
    return summary

if __name__ == "__main__":
    main()
