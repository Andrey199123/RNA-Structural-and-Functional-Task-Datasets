from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from .base import PipelineContext, TaskPipeline
from ..types import TaskState, TaskStatus
from ..validator import DataValidator


class Task4RBPBindingPipeline(TaskPipeline):
    task_id = "task4_rbp_binding"
    task_name = "Task 4: RNA-Protein Binding (RBP)"

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(ctx, output_dir=output_dir, config=config)
        self.task_dir = ctx.eval_root / "Task4"

    def run(self, status: TaskStatus) -> TaskStatus:
        script = self.task_dir / "prepare_task4_rbp_binding.py"
        if not script.exists():
            status.state = TaskState.FAILED
            status.message = f"Missing script: {script}"
            return status

        import subprocess
        import sys

        self.logger.info("Running %s", script)
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(self.task_dir),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "Task4 script failed")

        # Validation
        summary_tsv = self.task_dir / "task4_summary.tsv"
        h5_dir = self.task_dir / "hdf5"
        n_rows = 0
        if summary_tsv.exists():
            import csv

            with summary_tsv.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                n_rows = sum(1 for _ in reader)

        n_h5 = len(list(h5_dir.glob("*.h5"))) if h5_dir.exists() else 0
        v = DataValidator()
        checks = [
            v.validate_counts(observed=n_rows, expected=17, tolerance_ratio=0.0, name="rbp_rows"),
            v.validate_counts(observed=n_h5, expected=17, tolerance_ratio=0.0, name="h5_files"),
        ]
        report_path = self.task_dir / "validation.json"
        v.write_report(report=v.report(checks), path=report_path)

        status.output_dir = str(self.task_dir)
        status.artifacts = [
            str(self.task_dir / "task4_summary.tsv"),
            str(self.task_dir / "task4_README.md"),
            str(self.task_dir / "hdf5"),
            str(report_path),
        ]
        status.message = "Prepared PrismNet HDF5 files"
        return status
