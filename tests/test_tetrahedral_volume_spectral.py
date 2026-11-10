from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from tetrahedral_volume_backend import (
    apply_tetra_volume_tensor,
    tetra_volume_diagnostics,
)


class TetrahedralVolumeSpectralTests(unittest.TestCase):
    def test_raw_q_hermiticity_and_m_independence_are_measured_before_cleanup(self):
        result = tetra_volume_diagnostics((1, 1, 1, 1), 1)
        self.assertLess(max(result["raw_Q_hermiticity_defects"]), 1e-10)
        self.assertLess(result["max_raw_Q_M_block_defect"], 1e-10)
        self.assertLess(result["V4_minus_Q2_relative_error"], 1e-9)
        self.assertTrue(result["zero_aware_functional_calculus"])

    def test_complete_magnetic_tensor_is_reconstructed(self):
        rng = np.random.default_rng(20261110)
        for spins in ((1, 1, 1, 1), (0, 1, 1, 1)):
            tensor = rng.normal(size=tuple(s + 1 for s in spins))
            tensor = tensor + 1j * rng.normal(size=tuple(s + 1 for s in spins))
            output = apply_tetra_volume_tensor(tensor, spins)
            self.assertEqual(output.shape, tensor.shape)
            self.assertTrue(np.isfinite(output.real).all())
            self.assertTrue(np.isfinite(output.imag).all())


if __name__ == "__main__":
    unittest.main()
