#!/usr/bin/env python3
import unittest
import sympy as sp
from bqg_finite_master_kernel_gate import build_master_constraint


class FiniteMasterKernelGateTests(unittest.TestCase):
    def test_exact_common_kernel_and_projector(self):
        C1, C2 = sp.Matrix([[1, 0, 0]]), sp.Matrix([[0, 1, 0]])
        result = build_master_constraint([C1, C2], sp.eye(2))
        self.assertEqual(result["physical_kernel_dimension"], 1)
        self.assertEqual(result["master_rank"], 2)
        self.assertEqual(result["projector_matrix"], sp.diag(0, 0, 1))
        self.assertTrue(all(result["checks"].values()))

    def test_nontrivial_two_dimensional_physical_sector(self):
        C = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0]])
        result = build_master_constraint([C], sp.eye(1))
        self.assertEqual(result["physical_kernel_dimension"], 2)
        self.assertEqual(result["physical_projective_dimension"], 1)
        self.assertEqual(result["physical_fs_real_dimension"], 2)

    def test_positive_label_metric_with_off_diagonal_terms(self):
        C1, C2 = sp.Matrix([[1, 0]]), sp.Matrix([[1, -1]])
        result = build_master_constraint([C1, C2], sp.Matrix([[2, 1], [1, 2]]))
        self.assertEqual(result["physical_kernel_dimension"], 0)
        self.assertTrue(all(result["checks"].values()))

    def test_rejects_indefinite_metric(self):
        with self.assertRaises(ValueError):
            build_master_constraint([sp.eye(2)], sp.diag(1, -1))

    def test_rejects_nonhermitian_metric(self):
        with self.assertRaises(ValueError):
            build_master_constraint([sp.eye(2)], sp.Matrix([[1, 1], [0, 1]]))

    def test_rejects_incompatible_output_fibers(self):
        with self.assertRaises(ValueError):
            build_master_constraint([sp.eye(2), sp.ones(1, 2)], sp.eye(2))

    def test_rejects_float_input(self):
        with self.assertRaises(ValueError):
            build_master_constraint([sp.Matrix([[0.5, 0]])], sp.eye(1))


if __name__ == "__main__":
    unittest.main()
