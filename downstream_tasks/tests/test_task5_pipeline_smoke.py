from __future__ import annotations

import pickle
import tempfile
import unittest
from pathlib import Path

import numpy as np

from downstream_tasks.pipelines.base import PipelineContext
from downstream_tasks.pipelines.task5_rri import Task5RRIPipeline
from downstream_tasks.types import TaskState


class TestTask5PipelineSmoke(unittest.TestCase):
    def test_task5_pipeline_runs_with_tiny_pickles(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            logs = eval_root / "logs"
            ctx = PipelineContext(repo_root=root, eval_root=eval_root, logs_dir=logs)

            data_dir = eval_root / "Task5" / "data" / "DeepMirTar_SdA" / "data"
            data_dir.mkdir(parents=True, exist_ok=True)

            pos = {"p1": {"f1": 1.0, "f2": 2.0}, "p2": {"f2": 1.0}}
            neg = {"n1": {"f1": 0.0}, "n2": {"f1": 0.5, "f3": 3.0}}

            (data_dir / "Mark_features_dic_all.pkl").write_bytes(pickle.dumps(pos, protocol=pickle.HIGHEST_PROTOCOL))
            (data_dir / "Mark_N_features_dic_all.pkl").write_bytes(pickle.dumps(neg, protocol=pickle.HIGHEST_PROTOCOL))

            p = Task5RRIPipeline(ctx, config={"resume": False, "seed": 1})
            status = p.execute()
            self.assertEqual(status.state, TaskState.SUCCESS)

            out_npz = eval_root / "Task5" / "preprocessed" / "deepmirtar_features.npz"
            self.assertTrue(out_npz.exists())
            with np.load(out_npz, allow_pickle=True) as z:
                self.assertIn("X", z)
                self.assertIn("y", z)


if __name__ == "__main__":
    unittest.main()

