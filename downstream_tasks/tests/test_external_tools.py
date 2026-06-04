from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from downstream_tasks.tools.cdhit_est import parse_clstr, recommended_word_size
from downstream_tasks.tools.dssr import extract_pairs


class TestExternalTools(unittest.TestCase):
    def test_recommended_word_size(self):
        self.assertEqual(recommended_word_size(0.95), 8)
        self.assertEqual(recommended_word_size(0.90), 8)
        self.assertEqual(recommended_word_size(0.88), 7)
        self.assertEqual(recommended_word_size(0.85), 6)
        self.assertEqual(recommended_word_size(0.80), 5)
        self.assertEqual(recommended_word_size(0.75), 4)

    def test_parse_clstr(self):
        content = (
            ">Cluster 0\n"
            "0 100nt, >seqA... *\n"
            "1 100nt, >seqB... at 99.0%\n"
            ">Cluster 1\n"
            "0 50nt, >seqC... *\n"
        )
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.clstr"
            p.write_text(content, encoding="utf-8")
            clusters = parse_clstr(p)
            self.assertEqual(clusters, [["seqA", "seqB"], ["seqC"]])

    def test_extract_pairs(self):
        d = {"pairs": [{"nt1": "A.G1", "nt2": "A.C72"}, {"nt1_id": "B.U2", "nt2_id": "B.A3"}]}
        self.assertEqual(extract_pairs(d), [("A.G1", "A.C72"), ("B.U2", "B.A3")])


if __name__ == "__main__":
    unittest.main()

