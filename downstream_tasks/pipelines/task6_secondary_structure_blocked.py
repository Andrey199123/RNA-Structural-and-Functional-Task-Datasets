from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from .base import PipelineContext, TaskPipeline
from ..types import TaskState, TaskStatus


class Task6SecondaryStructureBlockedPipeline(TaskPipeline):
    task_id = "task6_secondary_structure"
    task_name = "Task 6: Secondary Structure Prediction (BLOCKED)"

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(ctx, output_dir=output_dir, config=config)

    def run(self, status: TaskStatus) -> TaskStatus:
        status.state = TaskState.BLOCKED
        status.message = "ON HOLD: requires explicit approval before implementation"
        return status
