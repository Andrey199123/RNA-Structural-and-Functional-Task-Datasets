from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .base import PipelineContext, TaskPipeline
from ..preprocessor import PreprocessorEngine
from ..types import TaskState, TaskStatus
from ..utils.files import write_json
from ..validator import DataValidator


def _load_pickle(path: Path):
    import pickle

    with path.open("rb") as f:
        return pickle.load(f, encoding="latin1")


class Task5RRIPipeline(TaskPipeline):
    task_id = "task5_rri"
    task_name = "Task 5: RNA-RNA Interaction (DeepMirTar features)"

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(ctx, output_dir=output_dir, config=config)
        self.task_dir = ctx.eval_root / "Task5"
        self.data_dir = self.task_dir / "data" / "DeepMirTar_SdA" / "data"
        self.prep_dir = self.task_dir / "preprocessed"

    def run(self, status: TaskStatus) -> TaskStatus:
        pos_path = self.data_dir / "Mark_features_dic_all.pkl"
        neg_path = self.data_dir / "Mark_N_features_dic_all.pkl"
        if not pos_path.exists() or not neg_path.exists():
            status.state = TaskState.SKIPPED
            status.message = f"Missing DeepMirTar pickle files in {self.data_dir}"
            return status

        pos = _load_pickle(pos_path)
        neg = _load_pickle(neg_path)

        # Build a stable feature order.
        feature_names = sorted({k for v in list(pos.values()) + list(neg.values()) for k in v.keys()})
        feat_index = {k: i for i, k in enumerate(feature_names)}

        ids: List[str] = []
        y: List[int] = []
        x = np.zeros((len(pos) + len(neg), len(feature_names)), dtype=np.float32)

        row = 0
        for key, feats in pos.items():
            ids.append(str(key))
            y.append(1)
            for fk, fv in feats.items():
                x[row, feat_index[fk]] = float(fv)
            row += 1
        for key, feats in neg.items():
            ids.append(str(key))
            y.append(0)
            for fk, fv in feats.items():
                x[row, feat_index[fk]] = float(fv)
            row += 1

        engine = PreprocessorEngine()
        seed = int(self.config.get("seed", 42))
        train_ratio = float(self.config.get("train_ratio", 0.8))
        val_ratio = float(self.config.get("val_ratio", 0.1))
        test_ratio = float(self.config.get("test_ratio", 0.1))
        train_idx, val_idx, test_idx = engine.split_stratified(
            y, train_ratio=train_ratio, val_ratio=val_ratio, test_ratio=test_ratio, seed=seed
        )

        self.prep_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            self.prep_dir / "deepmirtar_features.npz",
            ids=np.asarray(ids, dtype=object),
            X=x,
            y=np.asarray(y, dtype=np.int8),
            feature_names=np.asarray(feature_names, dtype=object),
            train_idx=train_idx,
            val_idx=val_idx,
            test_idx=test_idx,
        )

        meta = {
            "source": "DeepMirTar_SdA (feature pickles)",
            "note": "This repo snapshot contains engineered features, not raw paired sequences.",
            "n_samples": int(x.shape[0]),
            "n_features": int(x.shape[1]),
            "n_pos": int(sum(y)),
            "n_neg": int(len(y) - sum(y)),
            "split": {"train": int(len(train_idx)), "val": int(len(val_idx)), "test": int(len(test_idx))},
            "seed": seed,
            "ratios": {"train": train_ratio, "val": val_ratio, "test": test_ratio},
        }
        write_json(self.prep_dir / "summary.json", meta)

        v = DataValidator()
        train_ids = [ids[i] for i in train_idx.tolist()]
        val_ids = [ids[i] for i in val_idx.tolist()]
        test_ids = [ids[i] for i in test_idx.tolist()]
        checks = [
            v.validate_class_balance(y, min_fraction=0.01),
            v.validate_no_leakage(train_ids=train_ids, val_ids=val_ids, test_ids=test_ids),
            v.validate_counts(observed=int(x.shape[0]), expected=940, tolerance_ratio=0.0, name="n_samples_expected"),
        ]
        report = v.report(checks)
        report_path = self.prep_dir / "validation.json"
        v.write_report(report=report, path=report_path)

        status.output_dir = str(self.prep_dir)
        status.artifacts = [
            str(self.prep_dir / "deepmirtar_features.npz"),
            str(self.prep_dir / "summary.json"),
            str(report_path),
        ]
        status.message = f"Prepared DeepMirTar features (n={x.shape[0]} d={x.shape[1]})"
        return status
