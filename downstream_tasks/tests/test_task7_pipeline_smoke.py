from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from downstream_tasks.pipelines.base import PipelineContext
from downstream_tasks.pipelines.task7_long_range import Task7LongRangePipeline
from downstream_tasks.types import TaskState


class TestTask7PipelineSmoke(unittest.TestCase):
    def test_task7_pipeline_runs_with_tiny_jsonl(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            logs = eval_root / "logs"
            ctx = PipelineContext(repo_root=root, eval_root=eval_root, logs_dir=logs)

            data_dir = eval_root / "Task7" / "data" / "bprna_spot0"
            data_dir.mkdir(parents=True, exist_ok=True)
            tr0 = data_dir / "TR0.jsonl"
            rows = [
                {"id": "s1", "sequence": "AU", "secondary_structure": "()"},
                {"id": "s2", "sequence": "AUGC", "secondary_structure": "([)]"},
            ]
            tr0.write_text("\n".join(json.dumps(r) for r in rows) + "\n", encoding="utf-8")

            p = Task7LongRangePipeline(ctx, config={"resume": False, "max_length": 10, "long_range_min_sep": 1})
            status = p.execute()
            self.assertEqual(status.state, TaskState.SUCCESS)

            out_npz = eval_root / "Task7" / "preprocessed" / "TR0_long_range_pairs.npz"
            self.assertTrue(out_npz.exists())
            with np.load(out_npz, allow_pickle=True) as z:
                self.assertIn("pairs", z)


if __name__ == "__main__":
    unittest.main()

