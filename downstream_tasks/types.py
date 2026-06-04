from __future__ import annotations

import enum
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


class TaskState(str, enum.Enum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"
    BLOCKED = "BLOCKED"


@dataclass
class TaskStatus:
    task_id: str
    state: TaskState
    message: str = ""
    started_at_epoch: Optional[int] = None
    finished_at_epoch: Optional[int] = None
    output_dir: Optional[str] = None
    artifacts: List[str] = field(default_factory=list)
    extra: Dict[str, Any] = field(default_factory=dict)

    @staticmethod
    def now_epoch() -> int:
        return int(time.time())

