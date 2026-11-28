#!/usr/bin/env python3
import unittest
import sympy as sp

from bqg_finite_master_kernel_gate import build_master_constraint
from bqg_w5_seed_gram_gate import (
    SEED_K2_VALUES, exact_seed_gram, exact_full_seed_gram, run_gate,
)


class W5SeedGramGateTests(unittest.TestCase):
    def test_five_seed_intertwiners_are_exactly_orthonormal(self):
        gram = exact_seed_gram()
        self.assertEqual(gram, sp.eye(len(SEED_K2_VALUES)))
        self.assertTrue(all(entry.is_Rational for entry in gram))

    def test_full_k5_seed_network_gram_factorizes_to_identity(self):
        full = exact_full_seed_gram()
        self.assertEqual(full, sp.eye(len(SEED_K2_VALUES)))
        self.assertEqual(full.rank(), 5)

    def test_seed_gram_plugs_into_metric_aware_projector_interface(self):
        # Synthetic constraint used only to test Q plumbing, not a W5 constraint.
        Q = exact_seed_gram()
        C = sp.Matrix([[1, 0, 0, 0, 0]])
        result = build_master_constraint([C], sp.eye(1), input_metric=Q)
        self.assertEqual(result["input_metric"], sp.eye(5))
        self.assertEqual(result["physical_kernel_dimension"], 4)
        self.assertTrue(all(result["checks"].values()))

    def test_report_declares_identity_input_metric(self):
        report = run_gate()
        self.assertTrue(report["passed"])
        self.assertEqual(report["input_metric_Q"],
                         [["1", "0", "0", "0", "0"],
                          ["0", "1", "0", "0", "0"],
                          ["0", "0", "1", "0", "0"],
                          ["0", "0", "0", "1", "0"],
                          ["0", "0", "0", "0", "1"]])


if __name__ == "__main__":
    unittest.main()
