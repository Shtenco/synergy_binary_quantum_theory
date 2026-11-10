import unittest

import numpy as np

from scripts.bqg_lorentzian_w5_joint_gram_gate import compute_joint_gram


def encode(matrix):
    return [[[float(z.real), float(z.imag)] for z in row] for row in matrix]


def payload(A, B):
    return {
        "input_basis_K": [0, 2, 4, 6, 8],
        "output_basis_id": "sha256:shared-output-basis",
        "A_output_basis_id": "sha256:shared-output-basis",
        "B_output_basis_id": "sha256:shared-output-basis",
        "A": encode(A),
        "B": encode(B),
        "normalization": {"c_E": [1.0, 0.0], "c_L": [1.0, 0.0]},
        "normalization_source": "test-only explicit convention",
    }


class TestW5JointGramGate(unittest.TestCase):
    def setUp(self):
        self.A = np.array(
            [[1.0, 0.0, 0.0, 0.0, 0.0],
             [0.0, 2.0, 0.0, 0.0, 0.0]],
            dtype=complex,
        )

    def test_rectangular_maps_include_structural_right_kernel(self):
        data = payload(self.A, np.zeros_like(self.A))
        result = compute_joint_gram(data)
        diag = result["diagnostics"]
        self.assertEqual(diag["Hsum_rank_at_tolerance"], 2)
        self.assertEqual(diag["Hsum_kernel_dimension_at_tolerance"], 3)
        self.assertEqual(diag["stacked_rank_at_tolerance"], 2)
        self.assertEqual(diag["stacked_kernel_dimension_at_tolerance"], 3)
        self.assertEqual(len(result["spectra"]["Hsum_all_five_right_singular_values"]), 5)

    def test_euclidean_lorentzian_cancellation_is_detected(self):
        data = payload(self.A, -self.A)
        result = compute_joint_gram(data)
        diag = result["diagnostics"]
        self.assertEqual(diag["Hsum_rank_at_tolerance"], 0)
        self.assertEqual(diag["Hsum_kernel_dimension_at_tolerance"], 5)
        self.assertEqual(diag["stacked_kernel_dimension_at_tolerance"], 3)

    def test_implicit_normalization_is_rejected(self):
        data = payload(self.A, np.zeros_like(self.A))
        del data["normalization"]
        with self.assertRaisesRegex(ValueError, "normalization.c_E"):
            compute_joint_gram(data)

    def test_normalization_without_provenance_is_rejected(self):
        data = payload(self.A, np.zeros_like(self.A))
        del data["normalization_source"]
        with self.assertRaisesRegex(ValueError, "normalization_source"):
            compute_joint_gram(data)

    def test_output_basis_mismatch_is_rejected(self):
        data = payload(self.A, np.zeros_like(self.A))
        data["B_output_basis_id"] = "different-basis"
        with self.assertRaisesRegex(ValueError, "same ordered output basis"):
            compute_joint_gram(data)

    def test_empty_output_carrier_is_rejected(self):
        A = np.zeros((0, 5), dtype=complex)
        data = payload(A, A)
        with self.assertRaisesRegex(ValueError, "at least one output row"):
            compute_joint_gram(data)

    def test_nonfinite_coefficients_are_rejected(self):
        data = payload(self.A, np.zeros_like(self.A))
        data["normalization"]["c_L"] = [float("nan"), 0.0]
        with self.assertRaisesRegex(ValueError, "must be finite"):
            compute_joint_gram(data)


if __name__ == "__main__":
    unittest.main()
