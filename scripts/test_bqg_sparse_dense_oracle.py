#!/usr/bin/env python3
"""Independent tiny dense-oracle checks for sparse complex linear algebra.

These tests use a deterministic synthetic finite-dimensional operator. They
validate sparse key aggregation, dense reconstruction, and the Hilbert-space
adjoint identity. They are NOT an execution or validation of the production
K5 Hamiltonian, volume operator, or Peter-Weyl column.
"""
from __future__ import annotations

import unittest

import numpy as np

from bqg_sparse_pruning import coherent_aggregate


def synthetic_operator_entries():
    """Return duplicate sparse matrix-element contributions (row, col, amp)."""
    return [
        (0, 0, 1.0 + 0.5j),
        (0, 1, -0.25 + 0.75j),
        (0, 1, 0.125 - 0.25j),
        (1, 0, 0.5 - 0.125j),
        (1, 2, 1.0e-14 + 2.0e-14j),
        (2, 1, -0.75 + 0.5j),
        (2, 2, 0.25 - 0.5j),
        (2, 2, -0.25 + 0.5j),  # exact coherent cancellation
        (3, 4, 0.375 + 0.125j),
        (4, 3, -0.625 + 0.25j),
        (4, 4, 1.25 + 0.0j),
    ]


def dense_oracle(entries, dimension):
    """Independent dense assembly; intentionally does not call sparse helper."""
    matrix = np.zeros((dimension, dimension), dtype=np.complex128)
    for row, col, amplitude in entries:
        matrix[row, col] += amplitude
    return matrix


def sparse_matrix(entries):
    return coherent_aggregate(
        [((row, col), amplitude) for row, col, amplitude in entries]
    )


def sparse_apply(operator, vector, dimension):
    contributions = []
    for (row, col), amplitude in operator.items():
        contributions.append((row, amplitude * vector[col]))
    output = coherent_aggregate(contributions)
    return np.array([output.get(row, 0j) for row in range(dimension)])


class SparseDenseOracleTests(unittest.TestCase):
    def setUp(self):
        self.dimension = 5
        self.entries = synthetic_operator_entries()
        self.sparse = sparse_matrix(self.entries)
        self.dense = dense_oracle(self.entries, self.dimension)

    def test_sparse_assembly_matches_independent_dense_oracle(self):
        reconstructed = np.zeros_like(self.dense)
        for (row, col), amplitude in self.sparse.items():
            reconstructed[row, col] = amplitude
        np.testing.assert_allclose(
            reconstructed, self.dense, rtol=2e-15, atol=1e-30
        )

    def test_sparse_action_matches_dense_matrix_vector_product(self):
        vector = np.array(
            [0.5 + 0.25j, -1.0 + 0.5j, 0.125 - 0.75j,
             1.25 + 0.0j, -0.25 + 0.375j],
            dtype=np.complex128,
        )
        actual = sparse_apply(self.sparse, vector, self.dimension)
        expected = self.dense @ vector
        np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-29)

    def test_complex_adjoint_identity_against_dense_conjugate_transpose(self):
        x = np.array(
            [0.25 + 0.5j, -0.75 + 0.125j, 1.0 - 0.25j,
             0.5 + 0.75j, -0.125 - 0.5j],
            dtype=np.complex128,
        )
        y = np.array(
            [-0.5 + 0.25j, 0.125 + 1.0j, 0.75 - 0.5j,
             -1.0 + 0.125j, 0.375 + 0.625j],
            dtype=np.complex128,
        )
        lhs = np.vdot(x, self.dense @ y)
        rhs = np.vdot(self.dense.conj().T @ x, y)
        scale = max(1.0, abs(lhs), abs(rhs))
        self.assertLessEqual(abs(lhs - rhs) / scale, 5e-15)

    def test_tiny_nonzero_survives_and_exact_cancellation_disappears(self):
        self.assertIn((1, 2), self.sparse)
        self.assertEqual(self.sparse[(1, 2)], 1.0e-14 + 2.0e-14j)
        self.assertNotIn((2, 2), self.sparse)


if __name__ == "__main__":
    unittest.main(verbosity=2)
