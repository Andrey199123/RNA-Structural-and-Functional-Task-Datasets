from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple


RNA_ALPHABET = set("ACGU")


def normalize_rna(seq: str) -> str:
    return seq.strip().upper().replace("T", "U")


def is_valid_rna(seq: str, *, allow_n: bool = True) -> bool:
    seq = normalize_rna(seq)
    allowed = set("ACGUN") if allow_n else set("ACGU")
    return all(ch in allowed for ch in seq)


def filter_by_length(
    seqs: Iterable[str],
    *,
    exact_len: Optional[int] = None,
    min_len: Optional[int] = None,
    max_len: Optional[int] = None,
) -> List[str]:
    out: List[str] = []
    for s in seqs:
        s = normalize_rna(s)
        n = len(s)
        if exact_len is not None and n != exact_len:
            continue
        if min_len is not None and n < min_len:
            continue
        if max_len is not None and n > max_len:
            continue
        out.append(s)
    return out


@dataclass(frozen=True)
class Split:
    train_idx: Tuple[int, ...]
    val_idx: Tuple[int, ...]
    test_idx: Tuple[int, ...]

