from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .errors import ValidationError


@dataclass(frozen=True)
class ValidationCheck:
    name: str
    passed: bool
    details: str


@dataclass(frozen=True)
class ValidationReport:
    checks: Tuple[ValidationCheck, ...]

    @property
    def ok(self) -> bool:
        return all(c.passed for c in self.checks)


class DataValidator:
    def validate_counts(
        self,
        *,
        observed: int,
        expected: int,
        tolerance_ratio: float = 0.05,
        name: str = "sample_count",
    ) -> ValidationCheck:
        if expected <= 0:
            return ValidationCheck(name=name, passed=False, details="expected must be > 0")
        lo = int(round(expected * (1 - tolerance_ratio)))
        hi = int(round(expected * (1 + tolerance_ratio)))
        passed = lo <= observed <= hi
        return ValidationCheck(
            name=name,
            passed=passed,
            details=f"observed={observed} expected={expected} tolerance={tolerance_ratio:.2%} range=[{lo},{hi}]",
        )

    def validate_seq_lengths(
        self,
        seqs: Sequence[str],
        *,
        exact_len: Optional[int] = None,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
        name: str = "sequence_lengths",
    ) -> ValidationCheck:
        lengths = np.array([len(str(s)) for s in seqs], dtype=int)
        if len(lengths) == 0:
            return ValidationCheck(name=name, passed=False, details="no sequences provided")

        passed = True
        msgs: List[str] = [f"n={len(lengths)} min={lengths.min()} max={lengths.max()} mean={lengths.mean():.2f}"]
        if exact_len is not None:
            n_bad = int((lengths != exact_len).sum())
            if n_bad:
                passed = False
            msgs.append(f"exact_len={exact_len} bad={n_bad}")
        if min_len is not None:
            n_bad = int((lengths < min_len).sum())
            if n_bad:
                passed = False
            msgs.append(f"min_len={min_len} bad={n_bad}")
        if max_len is not None:
            n_bad = int((lengths > max_len).sum())
            if n_bad:
                passed = False
            msgs.append(f"max_len={max_len} bad={n_bad}")
        return ValidationCheck(name=name, passed=passed, details="; ".join(msgs))

    def validate_class_balance(
        self,
        labels: Sequence[int],
        *,
        min_fraction: float = 0.01,
        name: str = "class_distribution",
    ) -> ValidationCheck:
        y = np.asarray(labels, dtype=int)
        if y.size == 0:
            return ValidationCheck(name=name, passed=False, details="no labels provided")
        values, counts = np.unique(y, return_counts=True)
        total = counts.sum()
        fracs = {int(v): float(c) / float(total) for v, c in zip(values, counts)}
        passed = all(frac >= min_fraction for frac in fracs.values())
        return ValidationCheck(name=name, passed=passed, details=f"fractions={fracs}")

    def validate_no_leakage(
        self,
        *,
        train_ids: Iterable[str],
        val_ids: Iterable[str],
        test_ids: Iterable[str],
        name: str = "data_leakage",
    ) -> ValidationCheck:
        train_set = set(train_ids)
        val_set = set(val_ids)
        test_set = set(test_ids)
        inter_tv = train_set.intersection(val_set)
        inter_tt = train_set.intersection(test_set)
        inter_vt = val_set.intersection(test_set)
        passed = not (inter_tv or inter_tt or inter_vt)
        details = f"train∩val={len(inter_tv)} train∩test={len(inter_tt)} val∩test={len(inter_vt)}"
        return ValidationCheck(name=name, passed=passed, details=details)

    def write_report(self, *, report: ValidationReport, path) -> None:
        import json
        from pathlib import Path

        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "ok": report.ok,
            "checks": [{"name": c.name, "passed": c.passed, "details": c.details} for c in report.checks],
        }
        p.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def report(self, checks: Sequence[ValidationCheck]) -> ValidationReport:
        return ValidationReport(checks=tuple(checks))
