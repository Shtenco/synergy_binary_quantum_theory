import itertools
import unittest

import numpy as np


def random_operator_blocks(rng, dim, hermitian):
    raw = rng.normal(size=(2 * dim, 2 * dim)) + 1j * rng.normal(
        size=(2 * dim, 2 * dim)
    )
    op = raw + raw.conj().T if hermitian else raw - raw.conj().T
    return op.reshape(2, dim, 2, dim).transpose(0, 2, 1, 3)


def aux_trace_product(A, B, V):
    dim = A.shape[2]
    out = np.zeros((dim, dim), dtype=complex)
    for i, j, k in itertools.product(range(2), repeat=3):
        out += A[i, j] @ B[j, k] @ V[k, i]
    return out


def epsilon_triples(neighbors):
    out = []
    for r, omitted in enumerate(neighbors):
        tri = tuple(x for x in neighbors if x != omitted)
        for p in itertools.permutations(tri):
            positions = [tri.index(x) for x in p]
            inversions = sum(
                positions[i] > positions[j]
                for i in range(3) for j in range(i + 1, 3)
            )
            sign = (-1) ** r * (-1 if inversions % 2 else 1)
            out.append((sign, p))
    return out


class TestEpsilonKKVAdjoint(unittest.TestCase):
    def test_adjoint_is_reversed_order_partial_trace(self):
        rng = np.random.default_rng(20261110)
        dim = 3
        A = random_operator_blocks(rng, dim, False)
        B = random_operator_blocks(rng, dim, False)
        V = random_operator_blocks(rng, dim, True)

        lhs = aux_trace_product(A, B, V).conj().T
        rhs = aux_trace_product(V, B, A)
        np.testing.assert_allclose(lhs, rhs, rtol=2e-14, atol=2e-14)

    def test_operator_valued_auxiliary_trace_is_not_cyclic(self):
        rng = np.random.default_rng(20261111)
        dim = 3
        A = random_operator_blocks(rng, dim, False)
        B = random_operator_blocks(rng, dim, False)
        V = random_operator_blocks(rng, dim, True)

        lhs = aux_trace_product(A, B, V)
        cyclic_candidate = aux_trace_product(B, A, V)
        relative_difference = np.linalg.norm(lhs - cyclic_candidate) / max(
            np.linalg.norm(lhs), 1.0
        )
        self.assertGreater(relative_difference, 1e-3)

    def test_epsilon_orbit_has_24_unique_signed_triples(self):
        triples = epsilon_triples((1, 2, 3, 4))
        self.assertEqual(len(triples), 24)
        signs = {p: sign for sign, p in triples}
        self.assertEqual(len(signs), 24)
        for (a, b, c), sign in signs.items():
            self.assertEqual(signs[(b, a, c)], -sign)


if __name__ == "__main__":
    unittest.main()
