from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


def _load_task1_module():
    repo_root = Path(__file__).resolve().parents[2]
    path = repo_root / "downstream-tasks-evaluation" / "Task1" / "prepare_task1_rsa.py"
    spec = importlib.util.spec_from_file_location("prepare_task1_rsa", str(path))
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


class TestTask1RSATextParsing(unittest.TestCase):
    def test_parse_and_normalize(self):
        m = _load_task1_module()
        txt = "\n".join(
            [
                ">id1",
                "AUG",
                "0 200 400",
                ">id2",
                "AU",
                "10",
                "20",
            ]
        )
        recs = m.parse_rnasnap2_data_txt(txt)
        self.assertEqual([r["id"] for r in recs], ["id1", "id2"])
        self.assertEqual(recs[0]["sequence"], "AUG")
        self.assertEqual(recs[0]["asa_raw"], [0.0, 200.0, 400.0])
        # Per-base normalization: A/G=400, U/C=350 (clamped to 1.0)
        self.assertAlmostEqual(recs[0]["asa_normalized"][0], 0.0, places=7)
        self.assertAlmostEqual(recs[0]["asa_normalized"][1], 200.0 / 350.0, places=7)
        self.assertAlmostEqual(recs[0]["asa_normalized"][2], 1.0, places=7)
        self.assertEqual(recs[0]["label_mask"], [1, 1, 1])

    def test_parse_rnasol_benchmark_missing_mask(self):
        m = _load_task1_module()
        txt = "\n".join(
            [
                "# comment",
                ">id1",
                "AUG",
                "0 200 -1",
                "01-1",
            ]
        )
        recs = m.parse_rnasol_benchmark_txt(txt)
        self.assertEqual(len(recs), 1)
        self.assertEqual(recs[0]["id"], "id1")
        self.assertEqual(recs[0]["sequence"], "AUG")
        self.assertEqual(recs[0]["asa_raw"], [0.0, 200.0, -1.0])
        self.assertEqual(recs[0]["label_mask"], [1, 1, 0])
        self.assertAlmostEqual(recs[0]["asa_normalized"][0], 0.0, places=7)
        self.assertAlmostEqual(recs[0]["asa_normalized"][1], 200.0 / 350.0, places=7)
        self.assertAlmostEqual(recs[0]["asa_normalized"][2], 0.0, places=7)


if __name__ == "__main__":
    unittest.main()
