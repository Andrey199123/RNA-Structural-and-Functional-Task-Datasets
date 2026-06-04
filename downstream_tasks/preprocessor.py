from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np


class PreprocessorEngine:
    def filter_by_length(
        self,
        seqs: Sequence[str],
        *,
        exact_len: Optional[int] = None,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
    ) -> List[str]:
        out: List[str] = []
        for s in seqs:
            s = str(s).strip().upper().replace("T", "U")
            n = len(s)
            if exact_len is not None and n != exact_len:
                continue
            if min_len is not None and n < min_len:
                continue
            if max_len is not None and n > max_len:
                continue
            out.append(s)
        return out

    def split_random(
        self,
        n: int,
        *,
        train_ratio: float = 0.8,
        val_ratio: float = 0.1,
        test_ratio: float = 0.1,
        seed: int = 42,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if not np.isclose(train_ratio + val_ratio + test_ratio, 1.0):
            raise ValueError("Ratios must sum to 1.0")
        rng = np.random.default_rng(seed)
        idx = rng.permutation(n)
        n_train = int(round(train_ratio * n))
        n_val = int(round(val_ratio * n))
        train = idx[:n_train]
        val = idx[n_train : n_train + n_val]
        test = idx[n_train + n_val :]
        return train, val, test

    def split_stratified(
        self,
        labels: Sequence[int],
        *,
        train_ratio: float = 0.8,
        val_ratio: float = 0.1,
        test_ratio: float = 0.1,
        seed: int = 42,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        labels_arr = np.asarray(labels, dtype=int)
        classes = np.unique(labels_arr)
        rng = np.random.default_rng(seed)

        train_idx: List[int] = []
        val_idx: List[int] = []
        test_idx: List[int] = []

        for c in classes:
            idx = np.where(labels_arr == c)[0]
            idx = rng.permutation(idx)
            n = len(idx)
            n_train = int(round(train_ratio * n))
            n_val = int(round(val_ratio * n))
            train_idx.extend(idx[:n_train].tolist())
            val_idx.extend(idx[n_train : n_train + n_val].tolist())
            test_idx.extend(idx[n_train + n_val :].tolist())

        return (
            rng.permutation(np.array(train_idx, dtype=int)),
            rng.permutation(np.array(val_idx, dtype=int)),
            rng.permutation(np.array(test_idx, dtype=int)),
        )

    def one_hot(
        self,
        seqs: Sequence[str],
        *,
        alphabet: str = "ACGUN",
        pad_to: Optional[int] = None,
        dtype=np.float32,
    ) -> np.ndarray:
        mapping = {ch: i for i, ch in enumerate(alphabet)}
        k = len(alphabet)
        max_len = max(len(s) for s in seqs) if seqs else 0
        out_len = pad_to if pad_to is not None else max_len
        out = np.zeros((len(seqs), out_len, k), dtype=dtype)
        for row, s in enumerate(seqs):
            s = str(s).upper().replace("T", "U")
            for col, ch in enumerate(s[:out_len]):
                out[row, col, mapping.get(ch, mapping.get("N", k - 1))] = 1.0
        return out

    def zscore(self, x: np.ndarray, *, eps: float = 1e-8) -> np.ndarray:
        mean = x.mean(axis=0, keepdims=True)
        std = x.std(axis=0, keepdims=True)
        return (x - mean) / (std + eps)

    def minmax(self, x: np.ndarray, *, eps: float = 1e-8) -> np.ndarray:
        mn = x.min(axis=0, keepdims=True)
        mx = x.max(axis=0, keepdims=True)
        return (x - mn) / (mx - mn + eps)

