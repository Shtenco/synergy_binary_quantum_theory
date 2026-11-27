#!/usr/bin/env python3
import unittest
from bqg_w5_fixed_basis_operator import (
    assemble_columns, adjoint_columns, adjoint_residual,
    linearity_residuals, apply_columns,
)

class FixedBasisOperatorTests(unittest.TestCase):
    def setUp(self):
        self.basis = ("e0", "e1")
        self.matrix = {
            "e0": {"e0": 1+2j, "e1": 3-1j},
            "e1": {"e0": -2j, "e1": 0.5+0j},
        }
        _, self.cols, leakage = assemble_columns(self.basis, self.matrix.__getitem__)
        self.assertEqual(leakage, ())

    def test_complex_coefficients_and_tiny_nonzero_are_retained(self):
        tiny = 1e-200 + 2e-200j
        _, cols, leakage = assemble_columns(
            ("a",), lambda _: {"a": tiny, "outside": 0j}
        )
        self.assertEqual(leakage, (("0", "outside", 0j),) if False else ())
        self.assertEqual(cols[0][0], tiny)

    def test_fixed_basis_leakage_is_explicit(self):
        _, _, leakage = assemble_columns(("a",), lambda _: {"outside": 1e-40j})
        self.assertEqual(len(leakage), 1)
        self.assertEqual(leakage[0][1], "outside")
        self.assertEqual(leakage[0][2], 1e-40j)

    def test_linearity_of_matrix_action(self):
        report = linearity_residuals(
            self.cols, {0: 1+2j, 1: -0.25j}, {0: 2-1j, 1: 0.5},
            0.125-0.75j,
        )
        self.assertLess(report["homogeneity_absolute_residual"], 1e-14)
        self.assertLess(report["additivity_absolute_residual"], 1e-14)

    def test_direct_adjoint_matches_conjugate_transpose(self):
        direct = adjoint_columns(self.cols, 2)
        report = adjoint_residual(direct, self.cols, 2)
        self.assertEqual(report["absolute_frobenius_residual"], 0.0)

    def test_action(self):
        result = apply_columns(self.cols, {0: 1, 1: 2})
        self.assertAlmostEqual(result[0], 1+2j-4j)
        self.assertAlmostEqual(result[1], 3-1j+1)

if __name__ == "__main__":
    unittest.main()
