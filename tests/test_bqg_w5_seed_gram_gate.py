#!/usr/bin/env python3
import unittest
import sympy as sp

from bqg_finite_master_kernel_gate import build_master_constraint\nfrom bqg_w5_seed_gram_gate import (
    SEED_K2_VALUES, exact_seed_gram, run_gate,
)


class W5SeedGramGateTests(unittest.TestCase):
    def test_five_seed_intertwiners_are_exactly_orthonormal(self):
        gram = exact_seed_gram()
        self.assertEqual(gram, sp.eye(len(SEED_K2_VALUES)))
        self.assertTrue(all(entry.is_Rational for entry in gram))

    def test_seed_gram_plugs_into_metric_aware_projector_interface(self):\n        # Synthetic constraint used only to test Q plumbing, not a W5 constraint.\n        Q = exact_seed_gram()\n        C = sp.Matrix([[1, 0, 0, 0, 0]])\n        result = build_master_constraint([C], sp.eye(1), input_metric=Q)\n        self.assertEqual(result[\"input_metric\"], sp.eye(5))\n        self.assertEqual(result[\"physical_kernel_dimension\"], 4)\n        self.assertTrue(all(result[\"checks\"].values()))\n\n    def test_report_declares_identity_input_metric(self):
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
