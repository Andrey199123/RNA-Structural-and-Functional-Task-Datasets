from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from downstream_tasks.pipelines.base import PipelineContext
from downstream_tasks.pipelines.task3_contact_map import Task3ContactMapPipeline
from downstream_tasks.types import TaskState


_MINI_PDB = """\
ATOM      1  P     A A   1       0.000   0.000   0.000  1.00  0.00           P
ATOM      2  C4'   A A   1       0.000   0.000   1.500  1.00  0.00           C
ATOM      3  P     U A   2       0.000   0.000   5.000  1.00  0.00           P
ATOM      4  C4'   U A   2       0.000   0.000   6.500  1.00  0.00           C
TER
END
"""


class TestTask3PipelineSmoke(unittest.TestCase):
    def test_task3_pipeline_runs_without_network(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            logs = eval_root / "logs"
            ctx = PipelineContext(repo_root=root, eval_root=eval_root, logs_dir=logs)

            task_dir = eval_root / "Task3"
            data_dir = task_dir / "data"
            data_dir.mkdir(parents=True, exist_ok=True)

            # Prevent downloads/extraction by pre-creating expected paths.
            (data_dir / "train_set.tar.bz2").write_bytes(b"")
            (data_dir / "test_set.tar.bz2").write_bytes(b"")
            train_root = data_dir / "train_set" / "train_set"
            test_root = data_dir / "test_set" / "test_set"
            (train_root / "train_fasta").mkdir(parents=True, exist_ok=True)
            (train_root / "train_pdb").mkdir(parents=True, exist_ok=True)
            (test_root / "test_fasta").mkdir(parents=True, exist_ok=True)
            (test_root / "test_pdb").mkdir(parents=True, exist_ok=True)

            # One train and one test example.
            (train_root / "train_fasta" / "x1.fasta").write_text(">x1\nAU\n", encoding="utf-8")
            (train_root / "train_pdb" / "x1.pdb").write_text(_MINI_PDB, encoding="utf-8")
            (test_root / "test_fasta" / "x2.fasta").write_text(">x2\nAU\n", encoding="utf-8")
            (test_root / "test_pdb" / "x2.pdb").write_text(_MINI_PDB, encoding="utf-8")

            p = Task3ContactMapPipeline(ctx, config={"resume": False, "train_ratio": 1.0, "long_range_min_sep": 0})
            status = p.execute()
            self.assertEqual(status.state, TaskState.SUCCESS)

            out = task_dir / "preprocessed" / "train" / "x1.npz"
            self.assertTrue(out.exists())
            with np.load(out, allow_pickle=True) as z:
                cm = z["contact_map"]
                self.assertEqual(cm.shape, (2, 2))


if __name__ == "__main__":
    unittest.main()

