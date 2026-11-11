#!/usr/bin/env python3
"""Small-sector tests of the production K5 raw volume triple product.

These tests call the actual spin-matrix and volume routines from the K5
producer. They inspect the unsymmetrized Q before the historical Hermitian
completion, then test the derived sqrt(|Q|) matrix. They do not test the full
K5 HDA column or establish a continuum/physical volume theorem.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

import k5_peter_weyl_safe_hda_column as k5


class K5RawVolumeOperatorTests(unittest.TestCase):
    CASES = ((1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2))

    def test_raw_triple_product_is_hermitian_before_symmetrization(self):
        for spins in self.CASES:
            with self.subTest(spins=spins):
                q = k5.volume123_raw_operator(*spins)
                scale = max(np.linalg.norm(q), np.finfo(float).tiny)
                defect = np.linalg.norm(q - q.conj().T) / scale
                self.assertLessEqual(defect, 5e-14)

    def test_volume_matrix_is_hermitian_positive_semidefinite(self):
        for spins in self.CASES:
            with self.subTest(spins=spins):
                volume = k5.volume123_matrix(*spins)
                np.testing.assert_allclose(
                    volume, volume.conj().T, rtol=2e-13, atol=2e-13
                )
                eigenvalues = np.linalg.eigvalsh(volume)
                self.assertGreaterEqual(float(eigenvalues.min()), -2e-12)

    def test_raw_defect_ledger_is_populated_from_unsymmetrized_operator(self):
        for spins in self.CASES:
            with self.subTest(spins=spins):
                k5.volume123_matrix(*spins)
                record = k5.VOLUME_RAW_DEFECTS[spins]
                q = k5.volume123_raw_operator(*spins)
                expected = np.linalg.norm(q - q.conj().T) / max(
                    np.linalg.norm(q), np.finfo(float).tiny
                )
                self.assertAlmostEqual(
                    record["relative_antihermitian_defect"], expected, places=14
                )
                self.assertAlmostEqual(record["raw_norm"], np.linalg.norm(q), places=12)


    def test_raw_volume_is_covariant_under_all_leg_permutations(self):
        spins = (1, 2, 2)
        dims = tuple(s + 1 for s in spins)
        q = k5.volume123_raw_operator(*spins)
        old_indices = list(itertools.product(*(range(d) for d in dims)))
        old_flat = {index: n for n, index in enumerate(old_indices)}

        for permutation in itertools.permutations(range(3)):
            with self.subTest(permutation=permutation):
                new_spins = tuple(spins[i] for i in permutation)
                new_dims = tuple(dims[i] for i in permutation)
                new_indices = list(
                    itertools.product(*(range(d) for d in new_dims))
                )
                new_flat = {index: n for n, index in enumerate(new_indices)}
                unitary = np.zeros((len(new_indices), len(old_indices)), dtype=complex)
                for old_index in old_indices:
                    new_index = tuple(old_index[i] for i in permutation)
                    unitary[new_flat[new_index], old_flat[old_index]] = 1.0

                inversions = sum(
                    permutation[i] > permutation[j]
                    for i in range(3) for j in range(i + 1, 3)
                )
                parity = -1.0 if inversions % 2 else 1.0
                q_permuted = k5.volume123_raw_operator(*new_spins)
                np.testing.assert_allclose(
                    q_permuted,
                    parity * unitary @ q @ unitary.conj().T,
                    rtol=3e-13,
                    atol=3e-13,
                )


    def test_derived_volume_is_covariant_under_all_leg_permutations(self):
        # Q changes sign under odd permutations, but sqrt(|Q|) must not:
        # | -U Q U† | = U |Q| U†. This checks the production spectral
        # reconstruction rather than inferring it from raw-Q covariance.
        spins = (1, 2, 2)
        dims = tuple(s + 1 for s in spins)
        volume = k5.volume123_matrix(*spins)
        old_indices = list(itertools.product(*(range(d) for d in dims)))
        old_flat = {index: n for n, index in enumerate(old_indices)}

        for permutation in itertools.permutations(range(3)):
            with self.subTest(permutation=permutation):
                new_spins = tuple(spins[i] for i in permutation)
                new_dims = tuple(dims[i] for i in permutation)
                new_indices = list(
                    itertools.product(*(range(d) for d in new_dims))
                )
                new_flat = {index: n for n, index in enumerate(new_indices)}
                unitary = np.zeros((len(new_indices), len(old_indices)), dtype=complex)
                for old_index in old_indices:
                    new_index = tuple(old_index[i] for i in permutation)
                    unitary[new_flat[new_index], old_flat[old_index]] = 1.0

                volume_permuted = k5.volume123_matrix(*new_spins)
                np.testing.assert_allclose(
                    volume_permuted,
                    unitary @ volume @ unitary.conj().T,
                    rtol=4e-13,
                    atol=4e-13,
                )


    def test_oriented_volume_tensor_is_self_adjoint_at_every_vertex(self):
        # The vertex-dependent epsilon_j conjugations must preserve
        # self-adjointness of the local sqrt(|Q|) action. Test the actual
        # production tensor routine in the Frobenius inner product.
        rng = np.random.default_rng(20261111)
        for vertex in k5.VERT:
            for triple in self.CASES:
                spins_local = tuple(triple) + (1,)
                shape = tuple(s + 1 for s in spins_local)
                left = (
                    rng.normal(size=shape) + 1j * rng.normal(size=shape)
                )
                right = (
                    rng.normal(size=shape) + 1j * rng.normal(size=shape)
                )
                v_left = k5.apply_volume_tensor_oriented(
                    left, spins_local, vertex
                )
                v_right = k5.apply_volume_tensor_oriented(
                    right, spins_local, vertex
                )
                lhs = np.vdot(left, v_right)
                rhs = np.vdot(v_left, right)
                with self.subTest(vertex=vertex, spins=spins_local):
                    np.testing.assert_allclose(
                        lhs, rhs, rtol=3e-12, atol=3e-12,
                        err_msg=(
                            "oriented production volume tensor is not "
                            "self-adjoint in the local tensor inner product"
                        ),
                    )

    def test_stable_volume_norm_handles_large_finite_components(self):
        value = 1.2e308 + 1.2e308j
        observed = k5.stable_complex_frobenius_norm(np.array([[value]]))
        expected = np.hypot(1.2e308, 1.2e308)
        self.assertAlmostEqual(observed / expected, 1.0, places=14)

    def test_stable_volume_norm_rejects_nonfinite_components(self):
        with self.assertRaisesRegex(FloatingPointError, "non-finite"):
            k5.stable_complex_frobenius_norm(
                np.array([[complex(float("inf"), 0.0)]])
            )

    def test_volume_spectrum_resolution_is_scale_invariant(self):
        base = np.array([-2.0, -1e-16, 0.0, 3.0])
        observed = [
            k5.volume_spectrum_resolution(base * scale)
            for scale in (1e-320, 1e-120, 1e-12, 1.0, 1e12, 1e120)
        ]
        self.assertEqual(
            [item["resolved_numerical_nullity"] for item in observed],
            [2] * 6,
        )
        gaps = [item["normalized_gap_to_zero"] for item in observed]
        np.testing.assert_allclose(gaps, np.full(6, 2.0 / 3.0), rtol=2e-15)

    def test_volume_spectrum_resolution_handles_zero_operator(self):
        result = k5.volume_spectrum_resolution(np.zeros(4))
        self.assertEqual(result["spectral_scale"], 0.0)
        self.assertEqual(result["resolved_numerical_nullity"], 4)
        self.assertIsNone(result["min_resolved_abs_eigenvalue"])
        self.assertIsNone(result["max_scalar_sqrt_abs_derivative_on_resolved_modes"])

    def test_symmetrization_operator_error_bound_is_recorded(self):
        spins = (1, 2, 2)
        k5.volume123_matrix(*spins)
        q = k5.volume123_raw_operator(*spins)
        record = k5.VOLUME_RAW_DEFECTS[spins]
        exact_frobenius_delta = 0.5 * k5.stable_complex_frobenius_norm(
            q - q.conj().T
        )
        self.assertAlmostEqual(
            record["symmetrization_delta_frobenius_norm"],
            exact_frobenius_delta,
            places=14,
        )
        self.assertGreaterEqual(
            record["symmetrization_delta_spectral_norm_upper_bound"],
            record["symmetrization_delta_frobenius_norm"],
        )
        self.assertAlmostEqual(
            record["symmetrization_relative_frobenius_bound"],
            exact_frobenius_delta / record["raw_norm"],
            places=14,
        )

    def test_volume_square_reconstructs_absolute_raw_spectrum(self):
        # Check V^2 = |Q| from an independently assembled eigenspectral
        # reference. Positivity and Hermiticity alone do not establish that
        # the square-root spectral map was reconstructed correctly.
        for spins in self.CASES:
            with self.subTest(spins=spins):
                q_raw = k5.volume123_raw_operator(*spins)
                q_hermitian = (q_raw + q_raw.conj().T) / 2
                eigenvalues, eigenvectors = np.linalg.eigh(q_hermitian)
                abs_q = (
                    eigenvectors * np.abs(eigenvalues)
                ) @ eigenvectors.conj().T
                volume = k5.volume123_matrix(*spins)
                np.testing.assert_allclose(
                    volume @ volume,
                    abs_q,
                    rtol=5e-12,
                    atol=5e-12,
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
