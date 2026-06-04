from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from downstream_tasks.orchestrator import TaskPipelineOrchestrator
from downstream_tasks.config import dump_yaml


class TestOrchestratorConfig(unittest.TestCase):
    def test_load_task_configs_defaults_enabled(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            (eval_root / "config").mkdir(parents=True, exist_ok=True)
            orch = TaskPipelineOrchestrator(repo_root=root, eval_root=eval_root)
            cfg = orch.load_task_configs()
            # At least one task present and enabled
            self.assertTrue(any(v.get("enabled", False) for v in cfg.values()))

    def test_load_task_configs_merges_preprocessing(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            cfg_dir = eval_root / "config"
            cfg_dir.mkdir(parents=True, exist_ok=True)
            dump_yaml(cfg_dir / "task_configs.yaml", {"tasks": {"task1_rsa": {"enabled": True}}})
            dump_yaml(cfg_dir / "preprocessing_configs.yaml", {"task_defaults": {"seed": 7}, "tasks": {"task1_rsa": {"seed": 9}}})

            orch = TaskPipelineOrchestrator(repo_root=root, eval_root=eval_root)
            cfg = orch.load_task_configs()
            self.assertEqual(cfg["task1_rsa"]["seed"], 9)


if __name__ == "__main__":
    unittest.main()

