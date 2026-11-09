import itertools
import unittest

import numpy as np


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


class TestEpsilonKKVAdjoint(unittest.TestCase):
    def test_complete_24_term_block_is_antihermitian(self):
        rng = np.random.default_rng(20261110)
        dim = 3
        neighbors = (1, 2, 3, 4)
        triples = epsilon_triples(neighbors)
        self.assertEqual(len(triples), 24)
        self.assertEqual(len({p for _, p in triples}), 24)

        K = {e: random_operator_blocks(rng, dim, False) for e in neighbors}
        V = {e: random_operator_blocks(rng, dim, True) for e in neighbors}

        B = np.zeros((dim, dim), dtype=complex)
        for sign, (a, b, c) in triples:
            B += sign * aux_trace_product(K[a], K[b], V[c])

        scale = max(np.linalg.norm(B), 1.0)
        antihermitian_residual = np.linalg.norm(B + B.conj().T) / scale
        self.assertLess(antihermitian_residual, 2e-14)

        H = 1j * B
        hermitian_residual = np.linalg.norm(H - H.conj().T) / scale
        self.assertLess(hermitian_residual, 2e-14)

    def test_epsilon_sign_flips_under_first_two_leg_exchange(self):
        triples = dict(epsilon_triples((1, 2, 3, 4)))
        for (a, b, c), sign in triples.items():
            self.assertEqual(triples[(b, a, c)], -sign)


if __name__ == "__main__":
    unittest.main()
