from __future__ import annotations

import unittest

from downstream_tasks.validator import DataValidator


class TestDataValidator(unittest.TestCase):
    def test_validate_counts(self):
        v = DataValidator()
        ok = v.validate_counts(observed=100, expected=100, tolerance_ratio=0.0)
        self.assertTrue(ok.passed)
        bad = v.validate_counts(observed=90, expected=100, tolerance_ratio=0.05)
        self.assertFalse(bad.passed)

    def test_validate_no_leakage(self):
        v = DataValidator()
        chk = v.validate_no_leakage(train_ids=["a", "b"], val_ids=["c"], test_ids=["d"])
        self.assertTrue(chk.passed)
        chk2 = v.validate_no_leakage(train_ids=["a"], val_ids=["a"], test_ids=[])
        self.assertFalse(chk2.passed)


if __name__ == "__main__":
    unittest.main()

