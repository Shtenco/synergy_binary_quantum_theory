import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import bqg_k5_actual_master_projector_gate as gate


class TestActualK5MasterProjectorAlgebra(unittest.TestCase):
    def setUp(self):
        self.C = np.asarray([
            [[1, 0], [0, 0]],
            [[0, 1], [0, 0]],
            [[0, 0], [1, 0]],
        ], dtype=np.complex128)

    def test_identity_metric_is_positive_and_selects_common_kernel(self):
        M = gate.master_from_constraints(self.C, np.eye(3))
        values = np.linalg.eigvalsh(M)
        self.assertGreaterEqual(float(values.min()), -1e-12)
        values, P, gap, rank, _ = gate.spectral_projector(M)
        self.assertEqual(rank, 0)  # this test family has no common kernel
        self.assertGreater(gap, 0.0)
        self.assertLess(np.linalg.norm(P @ P - P), 1e-12)

    def test_positive_constraint_metric_preserves_common_kernel(self):
        C = np.asarray([
            [[1, 0, 0], [0, 0, 0]],
            [[0, 1, 0], [0, 0, 0]],
            [[1, 1, 0], [0, 0, 0]],
        ], dtype=np.complex128)
        X = np.asarray([[1, 1j, 0], [0, 2, 1j], [1, 0, 2]], dtype=np.complex128)
        G = X.conj().T @ X + 0.5 * np.eye(3)
        M0 = gate.master_from_constraints(C, np.eye(3))
        MG = gate.master_from_constraints(C, G)
        _, P0, _, r0, _ = gate.spectral_projector(M0)
        _, PG, _, rG, _ = gate.spectral_projector(MG)
        self.assertEqual(r0, 1)
        self.assertEqual(rG, 1)
        self.assertLess(np.linalg.norm(P0 - PG), 1e-8)

    def test_nonhermitian_metric_is_rejected(self):
        G = np.asarray([[1, 1], [0, 1]], dtype=np.complex128)
        with self.assertRaisesRegex(ValueError, "Hermitian"):
            gate.master_from_constraints(self.C[:2], G)

    def test_indefinite_metric_is_rejected(self):
        G = np.diag([1.0, -0.1, 1.0])
        with self.assertRaisesRegex(ValueError, "positive definite"):
            gate.master_from_constraints(self.C, G)

    def test_actual_constraint_assembler_has_frozen_dimensions(self):
        C, specs, target, dims = gate.assemble_constraint_family()
        self.assertEqual(C.shape, (60, 140, 32))
        self.assertEqual(len(specs), 60)
        self.assertEqual(target.shape, (32,))
        self.assertEqual(dims["gauss_hilbert_dimension"], 140)
        self.assertEqual(dims["active_domain_dimension"], 32)


if __name__ == "__main__":
    unittest.main(verbosity=2)
