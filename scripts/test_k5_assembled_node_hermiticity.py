#!/usr/bin/env python3
"""Finite-sector Hermiticity oracle for the assembled production K5 node H.

The test materializes a complete projected matrix P H_v P on two explicit
orthonormal recoupling kets. It independently obtains every matrix column
from the production no-prune H action and checks the full 2x2 matrix against
its conjugate transpose. Leakage outside the selected four-ket subspace is
reported, not silently treated as zero. This is a projected finite-sector
check, not a full-Hilbert-space or continuum theorem.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

import k5_peter_weyl_safe_hda_column as k5


class K5AssembledNodeHermiticityTests(unittest.TestCase):
    def test_assembled_node_H_has_hermitian_complete_projected_block(self):
        vertex = 0
        # Every edge starts at doubled spin S=1. Each T monomial hits any
        # edge at most twice, so Jmax2=S+2=3 is a sufficient upper wall for
        # this one-action calculation from the declared input sector.
        spin_wall = 3
        spins = (1,) * len(k5.EDGES)
        # Two distinct allowed, orthonormal K5 recoupling basis kets.
        input_keys = (
            (spins, (0, 0, 0, 0, 0)),
            (spins, (0, 0, 0, 0, 2)),
        )
        for key in input_keys:
            for v, channel in enumerate(key[1]):
                self.assertIn(
                    channel,
                    k5.allowed_k2_t(*k5.local_spins(key[0], v)),
                    msg=f"invalid recoupling channel {channel} at vertex {v}",
                )

        row = {key: i for i, key in enumerate(input_keys)}
        matrix = np.zeros((len(input_keys), len(input_keys)), dtype=complex)
        leakage = []
        for column, input_key in enumerate(input_keys):
            # This is the actual assembled H_v = 1/2 sum_sigma eps_sigma
            # (T_sigma + T_sigma^dagger), not a test-side reconstruction.
            result = k5.apply_H_cached_state(
                {input_key: 1.0 + 0.0j},
                vertex,
                spin_wall,
                mode="no_prune",
            )
            self.assertTrue(result, "assembled node-H column unexpectedly vanished")
            for output_key, amplitude in result.items():
                self.assertTrue(
                    math.isfinite(amplitude.real) and math.isfinite(amplitude.imag),
                    msg=f"non-finite H amplitude at {output_key}",
                )
                if output_key in row:
                    matrix[row[output_key], column] = amplitude

            omitted = [
                amplitude for output_key, amplitude in result.items()
                if output_key not in row
            ]
            leakage.append({
                "count": len(omitted),
                "l2": math.sqrt(math.fsum(abs(value) ** 2 for value in omitted)),
                "max_abs": max((abs(value) for value in omitted), default=0.0),
            })

        self.assertEqual(matrix.shape, (2, 2))
        defect = float(np.linalg.norm(matrix - matrix.conj().T))
        scale = max(float(np.linalg.norm(matrix)), np.finfo(float).tiny)
        relative_defect = defect / scale
        np.testing.assert_allclose(
            matrix,
            matrix.conj().T,
            rtol=2e-10,
            atol=2e-11,
            err_msg=(
                "complete projected production node-H matrix is not Hermitian"
            ),
        )
        print(
            "ASSEMBLED_NODE_H_PROJECTED_BLOCK",
            {
                "vertex": vertex,
                "basis_dimension": len(input_keys),
                "spin_wall_Jmax2": spin_wall,
                "matrix_norm": float(np.linalg.norm(matrix)),
                "absolute_hermiticity_defect": defect,
                "relative_hermiticity_defect": relative_defect,
                "leakage_outside_projected_basis": leakage,
            },
        )
        self.assertEqual(len(leakage), len(input_keys))
        for item in leakage:
            self.assertGreaterEqual(item["count"], 0)
            self.assertGreaterEqual(item["l2"], 0.0)
            self.assertGreaterEqual(item["max_abs"], 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
