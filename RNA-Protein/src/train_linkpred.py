#!/usr/bin/env python3
import os, sys, math, json, random
from pathlib import Path
import argparse
import yaml
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
from utils.io import load_npz_embeddings, pairs_to_index_tensors, set_seed

from sklearn.metrics import (
    matthews_corrcoef, accuracy_score, roc_auc_score, average_precision_score
)


def device_from_cfg(cfg):
    dev = cfg.get("runtime", {}).get("device", "auto")
    if dev == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(dev)


def build_edge_index(train_df, rna_ids, prot_ids):
    rmap = {k: i for i, k in enumerate(rna_ids)}
    pmap = {k: i for i, k in enumerate(prot_ids)}
    src, dst = [], []
    offset = len(rna_ids)
    for rna, prot in zip(train_df['rna_id'], train_df['prot_id']):
        if rna in rmap and prot in pmap:
            u = rmap[rna]; v = offset + pmap[prot]
            src.append(u); dst.append(v)
            src.append(v); dst.append(u)
    edge_index = torch.tensor([src, dst], dtype=torch.long)
    return edge_index


def sample_negatives(num_rna, num_prot, num_samples, existing_set, seed=42):
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
        out_u.append(u); out_v.append(v)
    return torch.tensor(out_u, dtype=torch.long), torch.tensor(out_v, dtype=torch.long)


def build_existing_edge_set(df, rna_ids, prot_ids):
    rmap = {k: i for i, k in enumerate(rna_ids)}
    pmap = {k: i for i, k in enumerate(prot_ids)}
    offset = len(rna_ids)
    s = set()
    for rna, prot in zip(df['rna_id'], df['prot_id']):
        if rna in rmap and prot in pmap:
            s.add((rmap[rna], offset + pmap[prot]))
    return s


def make_batches(pos_u, pos_v, batch_size):
    n = pos_u.shape[0]
    if n == 0:
        return
    idx = torch.randperm(n)
    for i in range(0, n, batch_size):
        j = min(i + batch_size, n)
        yield pos_u[idx[i:j]], pos_v[idx[i:j]]


def _metrics_from_arrays(y_true, y_prob, thr=0.5):
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    y_pred = (y_prob >= thr).astype(int)
    if y_true.size == 0:
        return {"acc": 0.0, "mcc": 0.0, "auroc": 0.0, "auprc": 0.0}
    return {
        "acc":  accuracy_score(y_true, y_pred),
        "mcc":  matthews_corrcoef(y_true, y_pred),
        "auroc": roc_auc_score(y_true, y_prob),
        "auprc": average_precision_score(y_true, y_prob),
    }


def _best_mcc_threshold(y_true, y_prob, grid=None):
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


def evaluate(model, head, x_rna, x_prot, edge_index, pos_df, rna_ids, prot_ids,
             eval_neg_ratio, existing_set, device):
    model.eval(); head.eval()
    with torch.no_grad():
        z = model(x_rna, x_prot, edge_index)  # [N_tot, D]
        num_rna = x_rna.shape[0]
        pu, pv = pairs_to_index_tensors(pos_df, rna_ids, prot_ids, device=device)
        pu = pu.cpu(); pv = pv.cpu()
        num_neg = int(math.ceil(len(pu) * eval_neg_ratio)) if len(pu) > 0 else 0
        neg_u, neg_v = sample_negatives(num_rna, x_prot.shape[0], num_neg, existing_set, seed=1234)
        pu = pu.to(device); pv = pv.to(device)
        neg_u = neg_u.to(device); neg_v = neg_v.to(device)
        pos_scores = torch.sigmoid(head(z[pu], z[pv])).cpu().numpy() if len(pu) > 0 else np.zeros((0,), dtype=np.float32)
        neg_scores = torch.sigmoid(head(z[neg_u], z[neg_v])).cpu().numpy() if len(neg_u) > 0 else np.zeros((0,), dtype=np.float32)

    y_prob = np.concatenate([pos_scores, neg_scores], axis=0)
    y_true = np.concatenate([np.ones_like(pos_scores), np.zeros_like(neg_scores)], axis=0)
    base = _metrics_from_arrays(y_true, y_prob, thr=0.5)
    base.update({"y_true": y_true, "y_prob": y_prob})
    return base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    cfg.setdefault("runtime", {})
    cfg.setdefault("train", {})
    cfg.setdefault("model", {})
    cfg["runtime"].setdefault("seed", 42)
    cfg["runtime"].setdefault("log_every", 1)
    cfg["runtime"].setdefault("ckpt_dir", "model_saves")
    cfg["runtime"].setdefault("device", "auto")
    cfg["train"].setdefault("epochs", 30)
    cfg["train"].setdefault("batch_size", 1024)
    cfg["train"].setdefault("lr", 1e-3)
    cfg["train"].setdefault("weight_decay", 1e-4)
    cfg["train"].setdefault("neg_ratio", 1.0)
    cfg["train"].setdefault("eval_neg_ratio", 1.0)
    cfg["train"].setdefault("early_stop_patience", 5)
    cfg["model"].setdefault("head_hidden", 256)

    set_seed(cfg["runtime"]["seed"])

    dd = cfg["data"]
    dd["pairs_train"] = dd["pairs_train"].format(**dd)
    dd["pairs_val"]   = dd["pairs_val"].format(**dd)
    dd["pairs_test"]  = dd["pairs_test"].format(**dd)

    device = device_from_cfg(cfg)
    ckpt_dir = Path(cfg["runtime"]["ckpt_dir"]); ckpt_dir.mkdir(parents=True, exist_ok=True)

    (rna_ids, Xr), (prot_ids, Xp) = load_npz_embeddings(cfg["embeddings"]["rna_npz"],
                                                        cfg["embeddings"]["prot_npz"])
    x_rna = Xr.to(device)
    x_prot = Xp.to(device)

    train_df = pd.read_csv(dd["pairs_train"]).drop_duplicates()
    val_df   = pd.read_csv(dd["pairs_val"]).drop_duplicates()
    test_df  = pd.read_csv(dd["pairs_test"]).drop_duplicates()

    all_df = pd.concat([train_df, val_df, test_df], ignore_index=True)
    exist_all = build_existing_edge_set(all_df, rna_ids, prot_ids)

    edge_index = build_edge_index(train_df, rna_ids, prot_ids).to(device)

    pos_u_tr, pos_v_tr = pairs_to_index_tensors(train_df, rna_ids, prot_ids, device=device)

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
    opt = optim.Adam(params, lr=cfg["train"]["lr"], weight_decay=cfg["train"]["weight_decay"])
    loss_fn = nn.BCEWithLogitsLoss()

    best_val_mcc = -1.0
    best_snapshot_path = Path(cfg["runtime"]["ckpt_dir"]) / f"{cfg.get('experiment_name','baseline')}.pt"
    patience = cfg["train"]["early_stop_patience"]
    bad_epochs = 0

    epochs = cfg["train"]["epochs"]
    bsz = cfg["train"]["batch_size"]
    neg_ratio = cfg["train"]["neg_ratio"]
    eval_neg_ratio = cfg["train"]["eval_neg_ratio"]

    for epoch in trange(1, epochs + 1, desc="epochs"):
        model.train(); head.train()
        total_loss = 0.0
        for pu, pv in make_batches(pos_u_tr, pos_v_tr, bsz):
            num_neg = int(math.ceil(len(pu) * neg_ratio))
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

        if epoch % cfg["runtime"]["log_every"] == 0:
            val_raw = evaluate(
                model, head, x_rna, x_prot, edge_index, val_df,
                rna_ids, prot_ids, eval_neg_ratio, exist_all, device
            )
            test_raw = evaluate(
                model, head, x_rna, x_prot, edge_index, test_df,
                rna_ids, prot_ids, eval_neg_ratio, exist_all, device
            )

            thr, best_mcc_val = _best_mcc_threshold(val_raw["y_true"], val_raw["y_prob"])
            val_m  = _metrics_from_arrays(val_raw["y_true"],  val_raw["y_prob"],  thr=thr)
            test_m = _metrics_from_arrays(test_raw["y_true"], test_raw["y_prob"], thr=thr)

            print(
                f"\nEpoch {epoch} | loss={total_loss:.4f} "
                f"| VAL acc={val_m['acc']:.4f} mcc={val_m['mcc']:.4f} auroc={val_m['auroc']:.4f} auprc={val_m['auprc']:.4f} "
                f"| TEST acc={test_m['acc']:.4f} mcc={test_m['mcc']:.4f} auroc={test_m['auroc']:.4f} auprc={test_m['auprc']:.4f} "
                f"| thr={thr:.2f}"
            )

            if val_m["mcc"] > best_val_mcc:
                best_val_mcc = val_m["mcc"]
                bad_epochs = 0
                torch.save({
                    "model": model.state_dict(),
                    "head": head.state_dict(),
                    "cfg": cfg,
                    "rna_ids": rna_ids,
                    "prot_ids": prot_ids,
                    "best_thr": float(thr),
                    "best_val_metrics": val_m,
                }, best_snapshot_path)
            else:
                bad_epochs += 1
                if bad_epochs >= patience:
                    print(f"Early stopping (no val MCC improvement for {patience} evals).")
                    break

    print(f"Best VAL MCC: {best_val_mcc:.4f} (checkpoint: {best_snapshot_path.name})")
    print("Training finished.")


if __name__ == "__main__":
    main()
