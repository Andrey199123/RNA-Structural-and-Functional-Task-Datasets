from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

from .base import PipelineContext, TaskPipeline
from ..types import TaskState, TaskStatus
from ..validator import DataValidator


class Task2SubcellularLocalizationPipeline(TaskPipeline):
    task_id = "task2_subcellular_localization"
    task_name = "Task 2: Subcellular Localization (DM3Loc/LOCAS; trunc8000)"

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(ctx, output_dir=output_dir, config=config)
        self.task_dir = ctx.eval_root / "Task2"

    def run(self, status: TaskStatus) -> TaskStatus:
        script = self.task_dir / "prepare_task2_subcellular_localization.py"
        if not script.exists():
            status.state = TaskState.FAILED
            status.message = f"Missing script: {script}"
            return status

        import subprocess
        import sys

        seed = int(self.config.get("seed", 42))
        folds = int(self.config.get("folds", 5))
        max_len = int(self.config.get("max_len", 8000))
        left = int(self.config.get("left", max_len // 2))
        right = int(self.config.get("right", max_len - left))
        keep_cytoplasm_only = bool(self.config.get("keep_cytoplasm_only", False))

        self.logger.info("Running %s", script)
        cmd = [
            sys.executable,
            str(script),
            "--seed",
            str(seed),
            "--folds",
            str(folds),
            "--max-len",
            str(max_len),
            "--left",
            str(left),
            "--right",
            str(right),
        ]
        if keep_cytoplasm_only:
            cmd.append("--keep-cytoplasm-only")
        result = subprocess.run(
            cmd,
            cwd=str(self.task_dir),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "Task2 script failed")

        prep = self.task_dir / "preprocessed"
        summary_path = prep / "summary.json"
        folds_path = prep / "folds.json"
        records_path = prep / "dm3loc_trunc8000.jsonl.gz"
        readme_path = prep / "README.md"

        if not summary_path.exists() or not folds_path.exists() or not records_path.exists() or not readme_path.exists():
            raise RuntimeError("Task2 preprocessed outputs missing; expected summary.json, folds.json, dm3loc_trunc8000.jsonl.gz, README.md")

        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        n = int(summary.get("counts", {}).get("n_sequences", 0))
        used_max = int(float(summary.get("lengths", {}).get("used", {}).get("max", 0)))

        v = DataValidator()
        expected_n = 17306 if keep_cytoplasm_only else 17277
        checks = [
            v.validate_counts(observed=n, expected=expected_n, tolerance_ratio=0.0, name="n_sequences_expected"),
            v.validate_counts(observed=used_max, expected=max_len, tolerance_ratio=0.0, name="max_len_used_expected"),
        ]
        report_path = prep / "validation.json"
        v.write_report(report=v.report(checks), path=report_path)

        status.output_dir = str(prep)
        status.artifacts = [
            str(records_path),
            str(folds_path),
            str(summary_path),
            str(readme_path),
            str(report_path),
        ]
        status.message = f"Prepared DM3Loc mRNA localization (n={n}, max_len={used_max})"
        return status
