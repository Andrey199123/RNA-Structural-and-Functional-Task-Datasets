from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from downstream_tasks.pipelines.base import PipelineContext, TaskPipeline
from downstream_tasks.types import TaskState, TaskStatus


class _DummyOk(TaskPipeline):
    task_id = "dummy_task"
    task_name = "Dummy"

    def run(self, status: TaskStatus) -> TaskStatus:  # type: ignore[override]
        status.message = "ran"
        return status


class _DummyFailIfRun(TaskPipeline):
    task_id = "dummy_task"
    task_name = "Dummy"

    def run(self, status: TaskStatus) -> TaskStatus:  # type: ignore[override]
        raise RuntimeError("should not run")


class TestPipelineResume(unittest.TestCase):
    def test_resume_skips_when_success(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            logs = eval_root / "logs"
            ctx = PipelineContext(repo_root=root, eval_root=eval_root, logs_dir=logs)

            first = _DummyOk(ctx, config={"resume": True})
            s1 = first.execute()
            self.assertEqual(s1.state, TaskState.SUCCESS)
            self.assertEqual(s1.message, "ran")

            second = _DummyFailIfRun(ctx, config={"resume": True})
            s2 = second.execute()
            self.assertEqual(s2.state, TaskState.SUCCESS)
            self.assertEqual(s2.message, "ran")

    def test_force_reruns_even_when_success(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            logs = eval_root / "logs"
            ctx = PipelineContext(repo_root=root, eval_root=eval_root, logs_dir=logs)

            first = _DummyOk(ctx, config={"resume": True})
            self.assertEqual(first.execute().state, TaskState.SUCCESS)

            second = _DummyFailIfRun(ctx, config={"resume": True, "force": True})
            s2 = second.execute()
            self.assertEqual(s2.state, TaskState.FAILED)


if __name__ == "__main__":
    unittest.main()

