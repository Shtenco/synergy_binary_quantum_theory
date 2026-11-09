import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_v10_w5_lorentzian_adjoint_closure as gate


def random_blocks(rng, dim, antihermitian):
    raw = rng.normal(size=(2 * dim, 2 * dim)) + 1j * rng.normal(
        size=(2 * dim, 2 * dim)
    )
    full = raw - raw.conj().T if antihermitian else raw + raw.conj().T
    return full.reshape(2, dim, 2, dim).transpose(0, 2, 1, 3)


def callback(blocks, antihermitian):
    def apply(state, source_v, edge, i, j, cutoff):
        vec = np.zeros(blocks[edge].shape[2], dtype=complex)
        for key, amp in state.items():
            vec[key[0]] += amp
        out_vec = blocks[edge][i, j] @ vec
        out = {
            (n, 0, 0): z
            for n, z in enumerate(out_vec)
            if abs(z) > 1e-14
        }
        if antihermitian:
            return out, {
                "outer_complete_basis_leakage": 0.0,
                "internal_volume_sector_leakage": 0.0,
                "complete_charge_basis_leakage": 0.0,
            }
        return out, 0.0
    return apply


def dense(state, dim):
    out = np.zeros(dim, dtype=complex)
    for key, amp in state.items():
        out[key[0]] += amp
    return out


class TestW5AdjointClosureAlgebra(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261110)
        self.dim = 4
        self.K = {e: random_blocks(rng, self.dim, True) for e in (1, 2, 3, 4)}
        self.V = {e: random_blocks(rng, self.dim, False) for e in (1, 2, 3, 4)}
        self.ck = callback(self.K, True)
        self.cv = callback(self.V, False)
        self.psi = {
            (0, 0, 0): 0.2 + 0.1j,
            (1, 0, 0): -0.3 + 0.2j,
            (2, 0, 0): 0.5 - 0.4j,
            (3, 0, 0): 0.1 + 0.3j,
        }

    def test_single_triple_direct_adjoint_matches_conjugate_transpose(self):
        forward, _, _ = gate.triple_forward(
            self.psi, 0, 1, 2, 3, 10, self.cv, self.ck
        )
        direct, _, _ = gate.triple_adjoint(
            self.psi, 0, 1, 2, 3, 10, self.cv, self.ck
        )
        x = dense(self.psi, self.dim)
        y = dense(forward, self.dim)
        z = dense(direct, self.dim)
        # Compare the two maps as matrices, not only on one vector.
        F = np.column_stack([
            dense(gate.triple_forward(
                {(n, 0, 0): 1.0}, 0, 1, 2, 3, 10, self.cv, self.ck
            )[0], self.dim)
            for n in range(self.dim)
        ])
        D = np.column_stack([
            dense(gate.triple_adjoint(
                {(n, 0, 0): 1.0}, 0, 1, 2, 3, 10, self.cv, self.ck
            )[0], self.dim)
            for n in range(self.dim)
        ])
        np.testing.assert_allclose(D, F.conj().T, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(z, F.conj().T @ x, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(y, F @ x, rtol=2e-13, atol=2e-13)

    def test_complete_epsilon_node_adjoint_matches_matrix_adjoint(self):
        triples = gate.epsilon_triples((1, 2, 3, 4))
        F = np.zeros((self.dim, self.dim), dtype=complex)
        D = np.zeros_like(F)
        for sign, (a, b, c) in triples:
            for n in range(self.dim):
                basis = {(n, 0, 0): 1.0}
                f, _, _ = gate.triple_forward(
                    basis, 0, a, b, c, 10, self.cv, self.ck
                )
                d, _, _ = gate.triple_adjoint(
                    basis, 0, a, b, c, 10, self.cv, self.ck
                )
                F[:, n] += sign * dense(f, self.dim)
                D[:, n] += sign * dense(d, self.dim)
        np.testing.assert_allclose(D, F.conj().T, rtol=5e-13, atol=5e-13)

    def test_epsilon_orbit_contains_24_unique_terms(self):
        triples = gate.epsilon_triples((1, 2, 3, 4))
        self.assertEqual(len(triples), 24)
        self.assertEqual(len({triple for _, triple in triples}), 24)

    def test_sparse_closed_matrix_adjoint_indexing(self):
        keys = [("basis", 0), ("basis", 1), ("basis", 2)]
        # Columns of a complex 3x3 forward matrix, represented sparsely.
        F = [
            {keys[0]: 1 + 2j, keys[2]: -0.5j},
            {keys[1]: 3 - 1j},
            {keys[0]: 2 + 0.25j, keys[1]: -4j},
        ]
        index = {key: i for i, key in enumerate(keys)}
        expected = gate._matrix_adjoint_expected_columns(F, keys, index)
        dense_F = np.column_stack([
            np.array([col.get(k, 0j) for k in keys], dtype=complex)
            for col in F
        ])
        dense_expected = np.column_stack([
            np.array([col.get(k, 0j) for k in keys], dtype=complex)
            for col in expected
        ])
        np.testing.assert_allclose(dense_expected, dense_F.conj().T)
        err, scale = gate.sparse_matrix_residual(expected, expected)
        self.assertEqual(err, 0.0)
        self.assertGreater(scale, 0.0)


if __name__ == "__main__":
    unittest.main()
