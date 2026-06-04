from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from downstream_tasks.orchestrator import TaskPipelineOrchestrator
from downstream_tasks.types import TaskState


class TestOrchestratorIntegration(unittest.TestCase):
    def test_read_existing_statuses_and_consolidate(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            (eval_root / "reports").mkdir(parents=True, exist_ok=True)

            # Write a single task status.
            tdir = eval_root / "reports" / "task1_rsa"
            tdir.mkdir(parents=True, exist_ok=True)
            (tdir / "status.json").write_text(
                json.dumps({"task_id": "task1_rsa", "state": "SUCCESS", "message": "ok"}, indent=2),
                encoding="utf-8",
            )

            orch = TaskPipelineOrchestrator(repo_root=root, eval_root=eval_root)
            statuses = orch.read_existing_statuses()
            self.assertEqual(statuses["task1_rsa"].state, TaskState.SUCCESS)
            # A different task should be pending (no status file).
            self.assertEqual(statuses["task3_contact_map"].state, TaskState.PENDING)

            out = orch.write_consolidated_report(statuses)
            self.assertTrue(out.exists())
            data = json.loads(out.read_text(encoding="utf-8"))
            self.assertIn("task1_rsa", data)


if __name__ == "__main__":
    unittest.main()

