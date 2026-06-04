from __future__ import annotations

import unittest
import numpy as np

from downstream_tasks.preprocessor import PreprocessorEngine


class TestPreprocessorEngine(unittest.TestCase):
    def test_filter_by_length_exact(self):
        pe = PreprocessorEngine()
        seqs = ["ACG", "ACGU", "A"]
        self.assertEqual(pe.filter_by_length(seqs, exact_len=3), ["ACG"])

    def test_split_stratified_preserves_classes(self):
        pe = PreprocessorEngine()
        labels = [0] * 50 + [1] * 50
        train, val, test = pe.split_stratified(labels, seed=1)
        self.assertEqual(set(np.unique(np.array(labels)[train])), {0, 1})
        self.assertEqual(set(np.unique(np.array(labels)[val])), {0, 1})
        self.assertEqual(set(np.unique(np.array(labels)[test])), {0, 1})

    def test_one_hot_shape(self):
        pe = PreprocessorEngine()
        x = pe.one_hot(["ACGUN"], pad_to=6)
        self.assertEqual(x.shape, (1, 6, 5))


if __name__ == "__main__":
    unittest.main()

