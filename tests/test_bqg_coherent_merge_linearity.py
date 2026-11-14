import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import peter_weyl_covariant_composition_gate as composition


class TestCoherentMergeLinearity(unittest.TestCase):
    def merge(self, terms, *, scale=1.0, tol=1e-11):
        out = {}
        for value in terms:
            composition.add(out, {"state": value}, scale=scale, tol=tol)
        return out.get("state", 0j)

    def test_subthreshold_partial_sum_is_not_deleted(self):
        # Each input is below TOL, but their coherent sum is above it.
        result = self.merge([0.6e-11, 0.6e-11])
        self.assertEqual(result, 1.2e-11)

    def test_small_nonzero_remainder_survives_cancellation(self):
        # Both contributions are individually significant; their remainder
        # is below TOL and must remain available to subsequent operators.
        result = self.merge([1.2e-11, -1.0e-11])
        self.assertAlmostEqual(result, 0.2e-11, delta=2e-27)
        self.assertNotEqual(result, 0j)

    def test_scaling_commutes_with_coherent_merge(self):
        terms = [1.2e-11, -1.0e-11]
        alpha = 3.0
        lhs = self.merge(terms, scale=alpha)
        rhs = alpha * self.merge(terms)
        self.assertAlmostEqual(lhs.real, rhs.real, delta=1e-26)
        self.assertAlmostEqual(lhs.imag, rhs.imag, delta=1e-26)

    def test_exact_cancellation_is_removed_but_not_subthreshold_remainder(self):
        exact = {}
        composition.add(exact, {"state": 2.0e-11})
        composition.add(exact, {"state": -2.0e-11})
        self.assertNotIn("state", exact)

        residual = {}
        composition.add(residual, {"state": 2.0e-11})
        composition.add(residual, {"state": -1.9e-11})
        self.assertAlmostEqual(residual["state"].real, 0.1e-11, delta=3e-28)

    def test_complex_real_and_imaginary_parts_are_accumulated(self):
        result = self.merge([0.4e-11 + 0.7e-11j, 0.4e-11 - 0.2e-11j])
        self.assertAlmostEqual(result.real, 0.8e-11, delta=1e-27)
        self.assertAlmostEqual(result.imag, 0.5e-11, delta=1e-27)


    def test_basis_relabeling_preserves_subthreshold_nonzero_amplitude(self):
        spins = tuple(range(10))
        labels = (
            ("C", 2, 0, 0, 0),
            ("G", 0),
            ("G", 0),
            ("G", 0),
            ("G", 0),
        )
        source = {(spins, labels): 1.0e-14 - 2.0e-14j}
        converted = composition.complete_to_covariant(source, source_v=0)
        self.assertEqual(len(converted), 1)
        amplitude = next(iter(converted.values()))
        self.assertEqual(amplitude, 1.0e-14 - 2.0e-14j)

    def test_invalid_audit_tolerance_is_rejected(self):
        with self.assertRaises(ValueError):
            composition.add({}, {"state": 1.0}, tol=float("nan"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
