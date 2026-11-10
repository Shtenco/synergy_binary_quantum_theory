import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_full_channel_flow_matrix as FLOW


class FullChannelFlowGuardTests(unittest.TestCase):
    def test_accepts_isolated_degenerate_master_channels(self):
        gap = FLOW.validate_master_spectrum(
            [5.0, 5.75, 22.5], [1e-14, 2e-14, 3e-14], "test"
        )
        self.assertAlmostEqual(gap, 0.75)

    def test_accepts_single_master_channel_without_interchannel_gap(self):
        self.assertIsNone(
            FLOW.validate_master_spectrum([11.8], [0.0], "single-channel test")
        )

    def test_rejects_split_two_dimensional_master_channel(self):
        with self.assertRaisesRegex(RuntimeError, "pair splitting"):
            FLOW.validate_master_spectrum([5.0], [1e-6], "split test")

    def test_rejects_unisolated_adjacent_master_channels(self):
        with self.assertRaisesRegex(RuntimeError, "not spectrally isolated"):
            FLOW.validate_master_spectrum([5.0, 5.0 + 1e-12], [0.0, 0.0], "gap test")

    def test_rejects_nonfinite_master_spectrum(self):
        with self.assertRaisesRegex(RuntimeError, "nonfinite"):
            FLOW.validate_master_spectrum([5.0, float("nan")], [0.0, 0.0], "nan test")

    def test_induced_refinement_tensor_reports_gram_and_singular_data(self):
        T = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]], dtype=complex)
        result = FLOW.refinement_tensor_diagnostics(T, "unit coisometry test")
        self.assertEqual(result["shape"], [2, 4])
        self.assertEqual(result["numerical_rank"], 2)
        self.assertLess(result["row_gram_identity_defect"], 1e-12)
        self.assertLess(result["column_gram_idempotence_defect"], 1e-12)
        self.assertIn("no unit-isometry assumption", result["normalization_claim"])

    def test_refinement_spectrum_checks_actual_gram_eigenvalues(self):
        T = np.array([[0.8, 0.0, 0.0], [0.0, 0.4, 0.0]], dtype=complex)
        result = FLOW.refinement_tensor_diagnostics(T, "strict spectral test")
        self.assertEqual(result["numerical_rank"], 2)
        self.assertEqual(result["singular_values"], [0.8, 0.4])
        self.assertLess(result["row_gram_singular_spectrum_defect"], 1e-12)
        self.assertLess(result["column_gram_singular_spectrum_defect"], 1e-12)
        self.assertLess(result["frobenius_singular_spectrum_defect"], 1e-12)
        self.assertAlmostEqual(result["frobenius_norm"] ** 2, 0.8**2 + 0.4**2)

    def test_induced_refinement_tensor_rejects_nonfinite_entries(self):
        T = np.array([[1.0, np.nan]], dtype=complex)
        with self.assertRaisesRegex(RuntimeError, "non-finite"):
            FLOW.refinement_tensor_diagnostics(T, "NaN test")

    def test_accepts_exactly_one_generated_two_dimensional_copy(self):
        rank, spread = FLOW.validate_generated_copy_spectrum(
            [0.5266343608, 0.5266343608, 1e-16, 1e-17], "copy test"
        )
        self.assertEqual(rank, 2)
        self.assertLess(spread, 1e-12)

    def test_rejects_generated_image_with_extra_rank(self):
        with self.assertRaisesRegex(RuntimeError, "rank 3"):
            FLOW.validate_generated_copy_spectrum(
                [0.5, 0.5, 0.01, 0.0], "extra-rank test"
            )

    def test_rejects_generated_image_with_rank_one(self):
        with self.assertRaisesRegex(RuntimeError, "rank 1"):
            FLOW.validate_generated_copy_spectrum(
                [0.5, 1e-15, 0.0, 0.0], "rank-one test"
            )

    def test_rejects_split_generated_copy_pair(self):
        with self.assertRaisesRegex(RuntimeError, "pair splitting"):
            FLOW.validate_generated_copy_spectrum(
                [0.5, 0.49, 1e-16, 0.0], "split-generated test"
            )


if __name__ == "__main__":
    unittest.main()
