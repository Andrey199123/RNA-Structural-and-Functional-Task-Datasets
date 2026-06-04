from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from downstream_tasks.pipelines.base import PipelineContext
from downstream_tasks.pipelines.task4_rbp_binding import Task4RBPBindingPipeline
from downstream_tasks.types import TaskState


class TestTask4PipelineSmoke(unittest.TestCase):
    def test_task4_pipeline_runs_with_stub_script(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            logs = eval_root / "logs"
            ctx = PipelineContext(repo_root=root, eval_root=eval_root, logs_dir=logs)

            task_dir = eval_root / "Task4"
            task_dir.mkdir(parents=True, exist_ok=True)

            # Stub script that writes expected artifacts without requiring PrismNet.
            script = task_dir / "prepare_task4_rbp_binding.py"
            script.write_text(
                (
                    "from pathlib import Path\n"
                    "td = Path(__file__).parent\n"
                    "(td / 'hdf5').mkdir(parents=True, exist_ok=True)\n"
                    "with (td / 'task4_summary.tsv').open('w', encoding='utf-8') as f:\n"
                    "    f.write('rbp\\tcount\\n')\n"
                    "    for i in range(17):\n"
                    "        f.write(f'RBP{i}\\t0\\n')\n"
                    "(td / 'task4_README.md').write_text('ok\\n', encoding='utf-8')\n"
                    "for i in range(17):\n"
                    "    (td / 'hdf5' / f'RBP{i}.h5').write_bytes(b'')\n"
                ),
                encoding="utf-8",
            )

            p = Task4RBPBindingPipeline(ctx, config={"resume": False})
            status = p.execute()
            self.assertEqual(status.state, TaskState.SUCCESS)
            self.assertTrue((task_dir / "validation.json").exists())


if __name__ == "__main__":
    unittest.main()

