from __future__ import annotations
import unittest
import numpy as np


class HDACompressionIdentityTests(unittest.TestCase):
    def test_exact_compression_identity_and_leakage_bound(self):
        rng = np.random.default_rng(20261119)
        n, k = 9, 3
        X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        A = (X + X.conj().T) / 2
        Y = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        B = (Y + Y.conj().T) / 2
        Z = rng.normal(size=(n, k)) + 1j * rng.normal(size=(n, k))
        W, _ = np.linalg.qr(Z)
        P = W @ W.conj().T
        Q = np.eye(n) - P
        Aeff, Beff = W.conj().T @ A @ W, W.conj().T @ B @ W
        lhs = W.conj().T @ (A @ B - B @ A) @ W
        rhs = Aeff @ Beff - Beff @ Aeff + W.conj().T @ A @ Q @ B @ W - W.conj().T @ B @ Q @ A @ W
        self.assertLess(np.linalg.norm(lhs - rhs), 1e-10)
        eta_a, eta_b = np.linalg.norm(Q @ A @ W, 2), np.linalg.norm(Q @ B @ W, 2)
        defect = np.linalg.norm((Aeff @ Beff - Beff @ Aeff) - lhs, 2)
        self.assertLessEqual(defect, 2 * eta_a * eta_b + 1e-10)

    def test_invariant_carrier_inherits_commutator_exactly(self):
        A = np.zeros((6, 6), dtype=complex)
        B = np.zeros((6, 6), dtype=complex)
        A[:3, :3] = np.array([[1, 1j, 0], [-1j, 2, 1], [0, 1, 3]], complex)
        B[:3, :3] = np.array([[2, 0, 1j], [0, 1, 1], [-1j, 1, 4]], complex)
        A[3:, 3:] = np.diag([5.0, 6.0, 7.0])
        B[3:, 3:] = np.diag([8.0, 9.0, 10.0])
        W = np.eye(6, dtype=complex)[:, :3]
        Aeff, Beff = W.conj().T @ A @ W, W.conj().T @ B @ W
        lhs = W.conj().T @ (A @ B - B @ A) @ W
        self.assertGreater(np.linalg.norm(Aeff @ Beff - Beff @ Aeff), 1e-6)
        self.assertLess(np.linalg.norm(Aeff @ Beff - Beff @ Aeff - lhs), 1e-12)


if __name__ == "__main__":
    unittest.main()
