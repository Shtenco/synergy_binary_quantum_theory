#!/usr/bin/env python3
"""Regression tests for deterministic coherent sparse amplitude aggregation."""
from __future__ import annotations

import itertools
import unittest
import numpy as np

from bqg_sparse_pruning import (
    PruningPolicy, aggregate_mappings, apply_pruning, coherent_aggregate,
    telescoping_error_bound, verify_ledger_identity,
)


class SparsePruningTests(unittest.TestCase):
    def test_permutation_invariance(self):
        terms = [(("k", 0), 1e-14 + 2e-14j),
                 (("k", 0), -3e-14 + 1e-14j),
                 (("k", 1), 2.0 - 1.0j)]
        expected = coherent_aggregate(terms)
        for perm in (terms, list(reversed(terms)), [terms[i] for i in (1, 2, 0)]):
            self.assertEqual(coherent_aggregate(perm), expected)

    def test_coherent_sum_before_threshold(self):
        tau = 1.0
        terms = [("x", 0.6*tau), ("x", -0.6*tau), ("x", 0.7*tau)]
        coherent = coherent_aggregate(terms)
        kept, ledger = apply_pruning(
            coherent, PruningPolicy("aggregate_then_prune", 0.5, "test")
        )
        self.assertEqual(kept, {"x": 0.7})
        self.assertEqual(ledger.discarded_keys, 0)

    def test_constructive_small_terms_survive_after_sum(self):
        terms = [("x", 0.4), ("x", 0.4), ("x", 0.4)]
        coherent = coherent_aggregate(terms)
        kept, _ = apply_pruning(
            coherent, PruningPolicy("aggregate_then_prune", 1.0, "test")
        )
        self.assertEqual(kept, {"x": 1.2000000000000002})

    def test_dense_oracle(self):
        columns = [
            {0: 1+2j, 2: 3-1j},
            {1: -2j, 2: 0.5+0.25j},
            {0: -1j, 3: 2+0j},
        ]
        aggregated = aggregate_mappings(columns)
        dense = np.zeros(4, dtype=np.complex128)
        for mapping in columns:
            for key, value in mapping.items():
                dense[key] += value
        got = np.zeros(4, dtype=np.complex128)
        for key, value in aggregated.items():
            got[key] = value
        np.testing.assert_allclose(got, dense, rtol=0, atol=1e-15)

    def test_pruning_ledger_identity(self):
        full = coherent_aggregate([("a", 1+1j), ("b", 0.01), ("c", -0.02j)])
        pruned, entry = apply_pruning(
            full, PruningPolicy("aggregate_then_prune", 0.05, "unit-test")
        )
        self.assertTrue(verify_ledger_identity(full, pruned, entry))

    def test_no_prune_never_discards(self):
        full = coherent_aggregate([("a", 1e-300), ("b", 0j)])
        kept, entry = apply_pruning(full, PruningPolicy())
        self.assertEqual(kept, full)
        self.assertEqual(entry.discarded_keys, 0)

    def test_telescoping_bound(self):
        # A0 then A1; perturbation at stage 0 is amplified by ||A1||.
        A1 = np.array([[2.0, 0.0], [0.0, 0.5]])
        d0 = np.array([0.01, -0.02])
        d1 = np.array([0.03, 0.0])
        exact_minus_approx = A1 @ d0 + d1
        bound = telescoping_error_bound(
            [np.linalg.norm(d0), np.linalg.norm(d1)],
            [np.linalg.norm(A1), 1.0],
        )
        self.assertLessEqual(np.linalg.norm(exact_minus_approx), bound + 1e-15)

    def test_cannot_prune_before_coherent_aggregation(self):
        with self.assertRaises(ValueError):
            apply_pruning(
                {"x": 0.1},
                PruningPolicy("aggregate_then_prune", 0.2, "bad"),
                coherent_aggregation_complete=False,
            )


if __name__ == "__main__":
    unittest.main()
