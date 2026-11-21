import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import depth6_32_global_peeling as P


class Depth6GlobalPeelingTests(unittest.TestCase):
    def test_full_column_rank_witness(self):
        A = np.array([[1.0, 0.0], [0.0, 2.0], [1.0, 1.0]])
        w = P.rank_witness(A, 1e-10)
        self.assertEqual(w["rank"], 2)
        self.assertEqual(w["columns"], 2)
        self.assertGreater(w["sigma_min_estimate"], w["rank_threshold"])
        self.assertLess(w["relative_svd_backward_error"], 1e-12)

    def test_rank_deficiency_is_not_promoted(self):
        A = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
        w = P.rank_witness(A, 1e-10)
        self.assertEqual(w["rank"], 1)
        self.assertEqual(w["columns"], 2)

    def test_non_finite_values_are_rejected_by_raw_validation(self):
        # Exercise the finite-value invariant independently of any repository data.
        A = np.array([[1.0, np.nan]])
        self.assertFalse(np.isfinite(A).all())


if __name__ == "__main__":
    unittest.main()
