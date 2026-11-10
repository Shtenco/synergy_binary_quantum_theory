import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_two_step_channel_composition as TWO


class TwoStepRefinementCompositionTests(unittest.TestCase):
    def test_composition_inserts_intermediate_projector_and_keeps_ancillas(self):
        T1 = np.arange(12, dtype=float).reshape(3, 2, 2).astype(complex)
        T2 = (np.arange(24, dtype=float).reshape(4, 3, 2) / 7).astype(complex)
        P = np.diag([1.0, 1.0, 0.0]).astype(complex)
        result = TWO.compose_refinement_tensors(T1, T2, P)
        self.assertEqual(result.shape, (4, 2, 2, 2))
        expected = np.einsum("rmb,mn,nca->rcab", T2, P, T1)
        self.assertLess(np.linalg.norm(result - expected), 1e-12)
        # Removing the intermediate projection is a genuinely different map.
        unprojected = np.einsum("rmb,mca->rcab", T2, T1)
        self.assertGreater(np.linalg.norm(result - unprojected), 1e-6)

    def test_rejects_nonprojector_intermediate_operator(self):
        T1 = np.zeros((2, 2, 2), dtype=complex)
        T2 = np.zeros((3, 2, 2), dtype=complex)
        P = np.array([[1.0, 0.2], [0.0, 1.0]], dtype=complex)
        with self.assertRaisesRegex(ValueError, "orthogonal projector"):
            TWO.compose_refinement_tensors(T1, T2, P)

    def test_refinement_spectral_distance_uses_rank_and_normalized_mass(self):
        left = {"singular_values": [2.0, 1.0], "numerical_rank": 2}
        right = {"singular_values": [1.0, 1.0], "numerical_rank": 2}
        result = TWO.normalized_refinement_spectral_distance(left, right)
        self.assertEqual(result["left_rank"], 2)
        self.assertEqual(result["right_rank"], 2)
        self.assertAlmostEqual(result["normalized_squared_singular_mass_l1"], 0.6)


if __name__ == "__main__":
    unittest.main()
