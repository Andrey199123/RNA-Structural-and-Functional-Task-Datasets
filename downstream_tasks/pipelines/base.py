from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

from ..logging_utils import setup_logging
from ..types import TaskState, TaskStatus
from ..utils.files import write_json


@dataclass(frozen=True)
class PipelineContext:
    repo_root: Path
    eval_root: Path
    logs_dir: Path


class TaskPipeline:
    task_id: str
    task_name: str

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        self.ctx = ctx
        self.output_dir = output_dir
        self.config = config or {}
        self.logger = setup_logging(name=f"downstream_tasks.{self.task_id}", log_dir=ctx.logs_dir)

    def status_path(self) -> Path:
        d = self.output_dir if self.output_dir is not None else (self.ctx.eval_root / "reports" / self.task_id)
        d.mkdir(parents=True, exist_ok=True)
        return d / "status.json"

    def write_status(self, status: TaskStatus) -> None:
        write_json(self.status_path(), status.__dict__)

    def read_status(self) -> Optional[TaskStatus]:
        path = self.status_path()
        if not path.exists():
            return None
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return None
        try:
            state = TaskState(str(data.get("state", TaskState.PENDING)))
        except Exception:
            state = TaskState.PENDING
        return TaskStatus(
            task_id=str(data.get("task_id", self.task_id)),
            state=state,
            message=str(data.get("message", "")),
            started_at_epoch=data.get("started_at_epoch"),
            finished_at_epoch=data.get("finished_at_epoch"),
            output_dir=data.get("output_dir"),
            artifacts=list(data.get("artifacts", []) or []),
            extra=dict(data.get("extra", {}) or {}),
        )

    def execute(self) -> TaskStatus:
        resume = bool(self.config.get("resume", True))
        force = bool(self.config.get("force", False))
        if resume and not force:
            existing = self.read_status()
            if existing is not None and existing.state == TaskState.SUCCESS:
                self.logger.info("Skipping (already SUCCESS): %s", self.task_name)
                return existing

        status = TaskStatus(task_id=self.task_id, state=TaskState.RUNNING, started_at_epoch=TaskStatus.now_epoch())
        try:
            self.logger.info("Starting: %s", self.task_name)
            status = self.run(status)
            if status.state == TaskState.RUNNING:
                status.state = TaskState.SUCCESS
            status.finished_at_epoch = TaskStatus.now_epoch()
            self.write_status(status)
            self.logger.info("Finished: %s (%s)", self.task_name, status.state)
            return status
        except BaseException as e:
            status.state = TaskState.FAILED
            status.message = f"{type(e).__name__}: {e}"
            status.finished_at_epoch = TaskStatus.now_epoch()
            self.write_status(status)
            self.logger.exception("Failed: %s", self.task_name)
            return status

    def run(self, status: TaskStatus) -> TaskStatus:
        raise NotImplementedError
