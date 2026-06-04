from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

from ..errors import ExternalToolError


@dataclass(frozen=True)
class CommandResult:
    argv: List[str]
    returncode: int
    stdout: str
    stderr: str


def which(program: str) -> Optional[str]:
    return shutil.which(program)


def require_executable(program: str) -> str:
    p = which(program)
    if not p:
        raise ExternalToolError(f"Executable not found on PATH: {program}")
    return p


def run(
    argv: Sequence[str],
    *,
    cwd: Optional[Path] = None,
    timeout_s: Optional[int] = None,
) -> CommandResult:
    result = subprocess.run(
        list(argv),
        cwd=str(cwd) if cwd is not None else None,
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    return CommandResult(argv=list(argv), returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)

