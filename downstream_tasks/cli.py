from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

from .orchestrator import TaskPipelineOrchestrator, DEFAULT_TASKS
from .types import TaskState
from .tools.subprocess_utils import which


def _orch_from_args(args: argparse.Namespace) -> TaskPipelineOrchestrator:
    repo_root = Path(args.repo_root).resolve() if getattr(args, "repo_root", None) else None
    eval_root = Path(args.eval_root).resolve() if getattr(args, "eval_root", None) else None
    config_dir = Path(args.config_dir).resolve() if getattr(args, "config_dir", None) else None
    return TaskPipelineOrchestrator(repo_root=repo_root, eval_root=eval_root, config_dir=config_dir)


def cmd_run_all(args: argparse.Namespace) -> int:
    orch = _orch_from_args(args)
    results = orch.execute_all_tasks(
        parallel=bool(args.parallel),
        max_workers=args.max_workers,
        resume=None if args.resume is None else bool(args.resume),
        force=bool(args.force),
    )
    orch.write_consolidated_report(results)
    orch.write_reproducibility_report()
    if args.json:
        print(json.dumps({k: v.__dict__ for k, v in results.items()}, indent=2))
    return 0 if all(r.state != TaskState.FAILED for r in results.values()) else 1


def cmd_run_task(args: argparse.Namespace) -> int:
    orch = _orch_from_args(args)
    status = orch.execute_single_task(
        args.task_id,
        resume=None if args.resume is None else bool(args.resume),
        force=bool(args.force),
    )
    if args.json:
        print(json.dumps(status.__dict__, indent=2))
    return 0 if status.state == TaskState.SUCCESS else 1


def cmd_status(args: argparse.Namespace) -> int:
    orch = _orch_from_args(args)
    if args.task_id:
        path = orch.get_task_status_path(args.task_id)
        if not path.exists():
            print(f"no status yet: {path}")
            return 1
        print(path.read_text(encoding="utf-8"))
        return 0

    for task_id in sorted(DEFAULT_TASKS):
        path = orch.get_task_status_path(task_id)
        if path.exists():
            print(f"{task_id}: {path}")
        else:
            print(f"{task_id}: (no status yet)")
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    orch = _orch_from_args(args)
    consolidated = orch.write_consolidated_report(orch.read_existing_statuses())
    repro = orch.write_reproducibility_report()
    print(f"wrote: {consolidated}")
    print(f"wrote: {repro}")
    return 0


def cmd_consolidate(args: argparse.Namespace) -> int:
    orch = _orch_from_args(args)
    statuses = orch.read_existing_statuses()
    out_path = orch.write_consolidated_report(statuses)
    if args.json:
        print(json.dumps({k: v.__dict__ for k, v in statuses.items()}, indent=2))
    else:
        print(f"wrote: {out_path}")
    return 0


def cmd_doctor(args: argparse.Namespace) -> int:
    checks = {
        "cd-hit-est": which("cd-hit-est"),
        "x3dna-dssr": which("x3dna-dssr"),
    }
    for name, path in checks.items():
        print(f"{name}: {path or 'NOT FOUND'}")
    return 0 if all(v is not None for v in checks.values()) else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m downstream_tasks.cli", description="Downstream tasks orchestrator")
    p.add_argument("--repo-root", default=None, help="Repository root (defaults to autodetect)")
    p.add_argument("--eval-root", default=None, help="downstream-tasks-evaluation directory (defaults to <repo>/downstream-tasks-evaluation)")
    p.add_argument("--config-dir", default=None, help="Config directory containing task_configs.yaml etc (defaults to <eval-root>/config)")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_all = sub.add_parser("run-all", help="Run all enabled tasks")
    p_all.add_argument("--json", action="store_true", help="Print JSON status output")
    p_all.add_argument("--parallel", action="store_true", help="Run enabled tasks in parallel (threaded)")
    p_all.add_argument("--max-workers", type=int, default=None, help="Max parallel workers (default: ThreadPoolExecutor default)")
    p_all.add_argument("--force", action="store_true", help="Force re-run even if status is already SUCCESS")
    p_all.add_argument(
        "--resume",
        dest="resume",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Resume from existing SUCCESS status (default: enabled)",
    )
    p_all.set_defaults(func=cmd_run_all)

    p_one = sub.add_parser("run-task", help="Run a single task")
    p_one.add_argument("task_id", choices=sorted(DEFAULT_TASKS.keys()))
    p_one.add_argument("--json", action="store_true", help="Print JSON status output")
    p_one.add_argument("--force", action="store_true", help="Force re-run even if status is already SUCCESS")
    p_one.add_argument(
        "--resume",
        dest="resume",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Resume from existing SUCCESS status (default: enabled)",
    )
    p_one.set_defaults(func=cmd_run_task)

    p_stat = sub.add_parser("status", help="Show task status")
    p_stat.add_argument("task_id", nargs="?", default=None)
    p_stat.set_defaults(func=cmd_status)

    p_rep = sub.add_parser("report", help="Write reproducibility report")
    p_rep.set_defaults(func=cmd_report)

    p_con = sub.add_parser("consolidate", help="Write consolidated report from existing status files")
    p_con.add_argument("--json", action="store_true", help="Print JSON status output")
    p_con.set_defaults(func=cmd_consolidate)

    p_doc = sub.add_parser("doctor", help="Check external tool availability")
    p_doc.set_defaults(func=cmd_doctor)

    return p


def main(argv: Any = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
