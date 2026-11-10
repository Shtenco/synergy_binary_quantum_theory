import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_full_channel_flow_matrix as FLOW
import bqg_two_step_channel_composition as COMPOSE


class TwoStepRefinementSpectrumTests(unittest.TestCase):
    def test_identical_normalized_spectra_have_zero_distance(self):
        diag = {
            "singular_values": [0.8, 0.4, 0.2],
        }
        result = COMPOSE.normalized_refinement_spectral_distance(diag, diag)
        self.assertAlmostEqual(result["normalized_squared_singular_mass_l1"], 0.0)
        self.assertAlmostEqual(result["normalized_squared_singular_mass_hellinger"], 0.0)

    def test_composed_tensor_has_actual_endpoint_spectrum(self):
        rng = np.random.default_rng(20261110)
        T1 = rng.normal(size=(3, 2, 2)) + 1j * rng.normal(size=(3, 2, 2))
        T2 = rng.normal(size=(4, 3, 2)) + 1j * rng.normal(size=(4, 3, 2))
        composed = np.einsum("rmb,mca->rcab", T2, T1, optimize=True)
        diag = FLOW.refinement_tensor_diagnostics(composed, "synthetic two-step tensor")
        self.assertEqual(diag["shape"], [4, 8])
        self.assertEqual(len(diag["singular_values"]), 4)
        self.assertLess(diag["row_gram_singular_spectrum_defect"], 1e-8)
        self.assertLess(diag["column_gram_singular_spectrum_defect"], 1e-8)


if __name__ == "__main__":
    unittest.main()
