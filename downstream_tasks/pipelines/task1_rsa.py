from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from .base import PipelineContext, TaskPipeline
from ..types import TaskState, TaskStatus


class Task1RSAPipeline(TaskPipeline):
    task_id = "task1_rsa"
    task_name = "Task 1: RNA Solvent Accessibility (RSA)"

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(ctx, output_dir=output_dir, config=config)
        self.task_dir = ctx.eval_root / "Task1"

    def run(self, status: TaskStatus) -> TaskStatus:
        script = self.task_dir / "prepare_task1_rsa.py"
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
            raise RuntimeError(result.stderr.strip() or "Task1 script failed")

        prep = self.task_dir / "preprocessed"
        status.output_dir = str(prep)
        status.artifacts = [
            str(prep / "train.json"),
            str(prep / "val.json"),
            str(prep / "test_ts44.json"),
            str(prep / "test_cn48.json"),
            str(prep / "statistics.json"),
            str(prep / "README.md"),
        ]
        status.message = "Prepared RSA splits (TR89 / TS44 / CN48)"
        return status
