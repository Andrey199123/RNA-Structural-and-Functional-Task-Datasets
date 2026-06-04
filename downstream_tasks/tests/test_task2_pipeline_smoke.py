from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from downstream_tasks.pipelines.base import PipelineContext
from downstream_tasks.pipelines.task2_subcellular_localization import Task2SubcellularLocalizationPipeline
from downstream_tasks.types import TaskState


class TestTask2PipelineSmoke(unittest.TestCase):
    def test_task2_pipeline_runs_with_stub_script(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            eval_root = root / "downstream-tasks-evaluation"
            logs = eval_root / "logs"
            ctx = PipelineContext(repo_root=root, eval_root=eval_root, logs_dir=logs)

            task_dir = eval_root / "Task2"
            task_dir.mkdir(parents=True, exist_ok=True)

            script = task_dir / "prepare_task2_subcellular_localization.py"
            script.write_text(
                (
                    "import gzip, json\n"
                    "from pathlib import Path\n"
                    "td = Path(__file__).parent\n"
                    "out = td / 'preprocessed'\n"
                    "out.mkdir(parents=True, exist_ok=True)\n"
                    "(out / 'folds.json').write_text(json.dumps({'k_folds': 5}), encoding='utf-8')\n"
                    "with gzip.open(out / 'dm3loc_trunc8000.jsonl.gz', 'wt', encoding='utf-8') as f:\n"
                    "    f.write(json.dumps({'id':'stub','sequence':'A','labels_6':[1,0,0,0,0,0]})+'\\n')\n"
                    "(out / 'summary.json').write_text(json.dumps({'counts': {'n_sequences': 17277}, 'lengths': {'used': {'max': 8000}}}, indent=2), encoding='utf-8')\n"
                    "(out / 'README.md').write_text('ok\\n', encoding='utf-8')\n"
                ),
                encoding="utf-8",
            )

            p = Task2SubcellularLocalizationPipeline(ctx, config={"resume": False})
            status = p.execute()
            self.assertEqual(status.state, TaskState.SUCCESS)
            self.assertTrue((task_dir / "preprocessed" / "validation.json").exists())


if __name__ == "__main__":
    unittest.main()
