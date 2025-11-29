from sklearn.metrics import matthews_corrcoef, roc_auc_score, average_precision_score, accuracy_score
import numpy as np

def binary_metrics(y_true, y_prob, thr=0.5):
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    y_pred = (y_prob >= thr).astype(int)

    acc = accuracy_score(y_true, y_pred)
    mcc = matthews_corrcoef(y_true, y_pred) if len(np.unique(y_true)) > 1 else 0.0
    auroc = roc_auc_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else 0.5
    auprc = average_precision_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else 0.5
    return dict(acc=acc, mcc=mcc, auroc=auroc, auprc=auprc)

def best_mcc_threshold(y_true, y_prob, grid=None):
    if grid is None:
        grid = np.linspace(0.05, 0.95, 37)
    best_t, best_m = 0.5, -1.0
    for t in grid:
        m = matthews_corrcoef((y_true>0.5).astype(int), (y_prob>=t).astype(int))
        if m > best_m: best_m, best_t = m, t
    return best_t, best_m

def metrics_with_threshold(y_true, y_prob, thr):
    pred = (y_prob >= thr).astype(int)
    return {
        "acc": accuracy_score(y_true, pred),
        "mcc": matthews_corrcoef(y_true, pred),
        "auroc": roc_auc_score(y_true, y_prob),
        "auprc": average_precision_score(y_true, y_prob),
    }
