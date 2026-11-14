#!/usr/bin/env python3
"""L1 characterization tests for thresholded production accumulators.

The remaining expected-failure tests expose known nonlinearity in the legacy
K-leg add_state and permutation-order paths. The COMP.add subthreshold
additivity check now passes and is a hard acceptance test.
"""
from __future__ import annotations

import itertools
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_sparse_linear_accumulator import (
    SparseLinearAccumulator,
    apply_sparse_matrix_columns,
    sparsify_matrix_columns,
)
import peter_weyl_covariant_composition_gate as COMP
import peter_weyl_covariant_K_leg_gate as CK


class TestStrictSparseLinearAccumulator(unittest.TestCase):
    def test_subthreshold_terms_survive_coherent_accumulation(self):
        tol = 1e-12
        accumulator = SparseLinearAccumulator()
        accumulator.add({"state": 0.6 * tol})
        accumulator.add({"state": 0.6 * tol})
        result = accumulator.finish()
        self.assertIn("state", result)
        self.assertAlmostEqual(result["state"].real / tol, 1.2, places=14)

    def test_complex_homogeneity_and_additivity(self):
        x = {"a": 1.25 + 0.5j, "b": -0.25j}
        y = {"a": -0.5 + 0.75j, "c": 2.0}
        alpha = -0.3 + 0.8j

        combined = SparseLinearAccumulator()
        combined.add(x)
        combined.add(y)
        expected_add = {"a": x["a"] + y["a"], "b": x["b"], "c": y["c"]}
        self.assertEqual(combined.finish(), expected_add)

        scaled = SparseLinearAccumulator()
        scaled.add(x, alpha)
        self.assertEqual(
            scaled.finish(),
            {key: alpha * value for key, value in x.items()},
        )

    def test_accumulation_is_stable_under_permutation(self):
        values = [
            1.1e-12 + 0.2e-12j,
            -0.9e-12 + 0.4e-12j,
            0.9e-12 - 0.1e-12j,
            0.6e-12 + 0.3e-12j,
        ]
        results = []
        for permutation in itertools.permutations(values):
            acc = SparseLinearAccumulator()
            for value in permutation:
                acc.add({"k": value})
            results.append(acc.finish()["k"])
        for value in results[1:]:
            self.assertAlmostEqual(value.real, results[0].real, places=27)
            self.assertAlmostEqual(value.imag, results[0].imag, places=27)

    def test_fixed_matrix_mask_preserves_linearity(self):
        columns = [
            {"r": 1.0 + 0.0j, "s": 1e-10 + 0j},
            {"r": 0.5j, "s": 2.0 + 0j},
        ]
        fixed = sparsify_matrix_columns(columns, 1e-9)
        x = [0.25 + 0.5j, -2.0 + 0.75j]
        y = [1.5 - 0.25j, 0.5 + 1.0j]
        lhs = apply_sparse_matrix_columns(fixed, [a + b for a, b in zip(x, y)])
        fx = apply_sparse_matrix_columns(fixed, x)
        fy = apply_sparse_matrix_columns(fixed, y)
        keys = set(lhs) | set(fx) | set(fy)
        for key in keys:
            self.assertAlmostEqual(
                lhs.get(key, 0j), fx.get(key, 0j) + fy.get(key, 0j), places=14
            )

    def test_nonfinite_scaled_amplitude_is_rejected(self):
        accumulator = SparseLinearAccumulator()
        with self.assertRaises(ValueError):
            accumulator.add({"huge": 1e308}, 1e308)

    def test_nonfinite_matrix_element_is_rejected(self):
        with self.assertRaises(ValueError):
            sparsify_matrix_columns([{"bad": complex(float("nan"), 0)}], 0.0)


class TestLegacyProductionAccumulatorLinearity(unittest.TestCase):
    def test_COMP_add_is_additive_below_threshold(self):
        tol = COMP.TOL
        x = {"same": 0.6 * tol}
        y = {"same": 0.6 * tol}
        combined = {}
        COMP.add(combined, {"same": x["same"] + y["same"]}, tol=tol)
        separate = {}
        COMP.add(separate, x, tol=tol)
        COMP.add(separate, y, tol=tol)
        self.assertEqual(combined, separate)

    @unittest.expectedFailure
    def test_CK_add_state_is_additive_below_threshold(self):
        tol = 1e-11
        x = {"same": 0.6 * tol}
        y = {"same": 0.6 * tol}
        combined = {}
        CK.add_state(combined, {"same": x["same"] + y["same"]}, tol=tol)
        separate = {}
        CK.add_state(separate, x, tol=tol)
        CK.add_state(separate, y, tol=tol)
        self.assertEqual(combined, separate)

    @unittest.expectedFailure
    def test_COMP_add_is_permutation_invariant(self):
        tol = COMP.TOL

        def fold(values):
            result = {}
            for value in values:
                COMP.add(result, {"same": value}, tol=tol)
            return result

        a = fold([1.1 * tol, -0.9 * tol, 0.9 * tol])
        b = fold([1.1 * tol, 0.9 * tol, -0.9 * tol])
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main(verbosity=2)
