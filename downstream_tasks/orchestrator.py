from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Type

from .pipelines.base import PipelineContext, TaskPipeline
from .pipelines.task1_rsa import Task1RSAPipeline
from .pipelines.task3_contact_map import Task3ContactMapPipeline
from .pipelines.task2_subcellular_investigation import Task2SubcellularInvestigationPipeline
from .pipelines.task2_subcellular_localization import Task2SubcellularLocalizationPipeline
from .pipelines.task4_rbp_binding import Task4RBPBindingPipeline
from .pipelines.task5_rri import Task5RRIPipeline
from .pipelines.task6_secondary_structure_blocked import Task6SecondaryStructureBlockedPipeline
from .pipelines.task7_long_range import Task7LongRangePipeline
from .types import TaskStatus
from .config import load_yaml
from .types import TaskState


DEFAULT_TASKS: Dict[str, Type[TaskPipeline]] = {
    Task1RSAPipeline.task_id: Task1RSAPipeline,
    Task2SubcellularInvestigationPipeline.task_id: Task2SubcellularInvestigationPipeline,
    Task2SubcellularLocalizationPipeline.task_id: Task2SubcellularLocalizationPipeline,
    Task3ContactMapPipeline.task_id: Task3ContactMapPipeline,
    Task4RBPBindingPipeline.task_id: Task4RBPBindingPipeline,
    Task5RRIPipeline.task_id: Task5RRIPipeline,
    Task6SecondaryStructureBlockedPipeline.task_id: Task6SecondaryStructureBlockedPipeline,
    Task7LongRangePipeline.task_id: Task7LongRangePipeline,
}


class TaskPipelineOrchestrator:
    def __init__(
        self,
        *,
        repo_root: Optional[Path] = None,
        eval_root: Optional[Path] = None,
        config_dir: Optional[Path] = None,
    ):
        self.repo_root = repo_root or Path(__file__).resolve().parents[1]
        self.eval_root = eval_root or (self.repo_root / "downstream-tasks-evaluation")
        self.logs_dir = self.eval_root / "logs"
        self.config_dir = config_dir or (self.eval_root / "config")
        self.ctx = PipelineContext(repo_root=self.repo_root, eval_root=self.eval_root, logs_dir=self.logs_dir)

    def load_task_configs(self) -> Dict[str, Dict[str, Any]]:
        task_cfg_path = self.config_dir / "task_configs.yaml"
        pre_cfg_path = self.config_dir / "preprocessing_configs.yaml"

        base_cfg: Dict[str, Any] = {}
        if task_cfg_path.exists():
            base_cfg = load_yaml(task_cfg_path)

        tasks = base_cfg.get("tasks", {}) if isinstance(base_cfg, dict) else {}
        if not isinstance(tasks, dict):
            raise ValueError(f"Invalid task_configs.yaml: expected tasks mapping, got {type(tasks).__name__}")

        pre_cfg: Dict[str, Any] = {}
        if pre_cfg_path.exists():
            pre_cfg = load_yaml(pre_cfg_path)
        defaults = pre_cfg.get("task_defaults", {}) if isinstance(pre_cfg, dict) else {}
        per_task = pre_cfg.get("tasks", {}) if isinstance(pre_cfg, dict) else {}
        if defaults is None:
            defaults = {}
        if per_task is None:
            per_task = {}
        if not isinstance(defaults, dict) or not isinstance(per_task, dict):
            raise ValueError("Invalid preprocessing_configs.yaml: expected task_defaults/tasks mappings")

        out: Dict[str, Dict[str, Any]] = {}
        for task_id in DEFAULT_TASKS:
            cfg = tasks.get(task_id, {})
            if cfg is None:
                cfg = {}
            if not isinstance(cfg, dict):
                raise ValueError(f"Invalid task config for {task_id}: expected mapping")
            merged = {"enabled": True, **defaults, **cfg}
            if task_id in per_task and isinstance(per_task[task_id], dict):
                merged.update(per_task[task_id])
            out[task_id] = merged
        return out

    def execute_single_task(self, task_id: str, *, resume: Optional[bool] = None, force: Optional[bool] = None) -> TaskStatus:
        if task_id not in DEFAULT_TASKS:
            raise KeyError(f"Unknown task_id: {task_id}")
        cfg = self.load_task_configs().get(task_id, {"enabled": True})
        if not cfg.get("enabled", True):
            return TaskStatus(task_id=task_id, state=TaskState.SKIPPED, message="Disabled by config")
        cfg = dict(cfg)
        if resume is not None:
            cfg["resume"] = bool(resume)
        if force is not None:
            cfg["force"] = bool(force)
        pipeline = DEFAULT_TASKS[task_id](self.ctx, config=cfg)
        return pipeline.execute()

    def execute_all_tasks(
        self,
        *,
        parallel: bool = False,
        max_workers: Optional[int] = None,
        resume: Optional[bool] = None,
        force: Optional[bool] = None,
    ) -> Dict[str, TaskStatus]:
        cfgs = self.load_task_configs()
        out: Dict[str, TaskStatus] = {}
        items = []
        for task_id, cls in DEFAULT_TASKS.items():
            cfg = cfgs.get(task_id, {"enabled": True})
            if not cfg.get("enabled", True):
                out[task_id] = TaskStatus(task_id=task_id, state=TaskState.SKIPPED, message="Disabled by config")
                continue
            cfg = dict(cfg)
            if resume is not None:
                cfg["resume"] = bool(resume)
            if force is not None:
                cfg["force"] = bool(force)
            items.append((task_id, cls, cfg))

        if parallel and items:
            from concurrent.futures import ThreadPoolExecutor, as_completed

            with ThreadPoolExecutor(max_workers=max_workers) as ex:
                futs = {ex.submit(cls(self.ctx, config=cfg).execute): task_id for task_id, cls, cfg in items}
                for fut in as_completed(futs):
                    task_id = futs[fut]
                    out[task_id] = fut.result()
        else:
            for task_id, cls, cfg in items:
                out[task_id] = cls(self.ctx, config=cfg).execute()
        return out

    def get_task_status_path(self, task_id: str) -> Path:
        return self.eval_root / "reports" / task_id / "status.json"

    def read_existing_statuses(self) -> Dict[str, TaskStatus]:
        import json

        out: Dict[str, TaskStatus] = {}
        for task_id in DEFAULT_TASKS:
            path = self.get_task_status_path(task_id)
            if not path.exists():
                out[task_id] = TaskStatus(task_id=task_id, state=TaskState.PENDING, message="No status yet")
                continue
            data = json.loads(path.read_text(encoding="utf-8"))
            try:
                state = TaskState(str(data.get("state", TaskState.PENDING)))
            except Exception:
                state = TaskState.PENDING
            out[task_id] = TaskStatus(
                task_id=task_id,
                state=state,
                message=str(data.get("message", "")),
                started_at_epoch=data.get("started_at_epoch"),
                finished_at_epoch=data.get("finished_at_epoch"),
                output_dir=data.get("output_dir"),
                artifacts=list(data.get("artifacts", []) or []),
                extra=dict(data.get("extra", {}) or {}),
            )
        return out

    def write_consolidated_report(self, statuses: Dict[str, TaskStatus]) -> Path:
        report_dir = self.eval_root / "reports"
        report_dir.mkdir(parents=True, exist_ok=True)
        out_path = report_dir / "consolidated_status.json"
        import json

        out_path.write_text(json.dumps({k: v.__dict__ for k, v in statuses.items()}, indent=2), encoding="utf-8")
        return out_path

    def write_reproducibility_report(self) -> Path:
        import json
        import platform
        import sys
        import time

        def _ver(mod: str) -> Optional[str]:
            try:
                m = __import__(mod)
                return getattr(m, "__version__", None)
            except Exception:
                return None

        report_dir = self.eval_root / "reports"
        report_dir.mkdir(parents=True, exist_ok=True)
        out_path = report_dir / "reproducibility.json"
        data = {
            "generated_at_epoch": int(time.time()),
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "packages": {
                "numpy": _ver("numpy"),
                "pandas": _ver("pandas"),
                "h5py": _ver("h5py"),
                "yaml": _ver("yaml"),
                "requests": _ver("requests"),
                "Bio": _ver("Bio"),
                "scipy": _ver("scipy"),
                "datasets": _ver("datasets"),
            },
        }
        out_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        return out_path
