#!/usr/bin/env python3
"""Non-vacuous finite-sector Hermiticity oracle for the assembled K5 node H.

The second basis ket is selected deterministically from a nonzero production
H action on the first ket. The complete projected 2x2 matrix therefore
contains a measured nonzero transition, rather than a potentially empty
same-spin block. The reverse matrix element is computed independently by
applying production H to the selected output ket. Leakage is reported.
This remains a two-state finite-sector test, not a full-sector or continuum
theorem.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

import k5_peter_weyl_safe_hda_column as k5


class K5AssembledNodeHermiticityTests(unittest.TestCase):
    def test_assembled_node_H_has_nonvacuous_hermitian_projected_block(self):
        vertex = 0
        spins = (1,) * len(k5.EDGES)
        input_key = (spins, (0, 0, 0, 0, 0))

        def validate_key(key):
            for v, channel in enumerate(key[1]):
                self.assertIn(
                    channel,
                    k5.allowed_k2_t(*k5.local_spins(key[0], v)),
                    msg=f"invalid recoupling channel {channel} at vertex {v}",
                )

        validate_key(input_key)

        # First H action: S=1, hence Jmax2=3 is above the proven S+2
        # reachable-spin ceiling for one action.
        forward_wall = max(input_key[0]) + 2
        forward = k5.apply_H_cached_state(
            {input_key: 1.0 + 0.0j},
            vertex,
            forward_wall,
            mode="no_prune",
        )
        self.assertTrue(forward, "assembled node-H column unexpectedly vanished")
        for output_key, amplitude in forward.items():
            self.assertTrue(
                math.isfinite(amplitude.real) and math.isfinite(amplitude.imag),
                msg=f"non-finite forward H amplitude at {output_key}",
            )

        candidates = [
            (key, amplitude)
            for key, amplitude in forward.items()
            if key != input_key and amplitude != 0j
        ]
        self.assertTrue(
            candidates,
            "H action produced no off-diagonal transition; projected test would be vacuous",
        )
        # Stable tie-break: largest |amplitude|, then lexicographic basis key.
        output_key, selected_amplitude = sorted(
            candidates, key=lambda item: (-abs(item[1]), item[0])
        )[0]
        validate_key(output_key)
        self.assertGreater(abs(selected_amplitude), 0.0)

        basis = (input_key, output_key)
        row = {key: i for i, key in enumerate(basis)}
        matrix = np.zeros((2, 2), dtype=complex)
        leakage = []

        # Column 0 is already the complete no-prune H(input_key) result.
        for output, amplitude in forward.items():
            if output in row:
                matrix[row[output], 0] = amplitude
        omitted = [
            amplitude for output, amplitude in forward.items() if output not in row
        ]
        leakage.append({
            "count": len(omitted),
            "l2": math.sqrt(math.fsum(abs(value) ** 2 for value in omitted)),
            "max_abs": max((abs(value) for value in omitted), default=0.0),
        })

        # Second action: the selected output has maximum doubled spin S_out;
        # Jmax2=S_out+2 avoids upper-wall clipping for this one H action.
        reverse_wall = max(output_key[0]) + 2
        reverse = k5.apply_H_cached_state(
            {output_key: 1.0 + 0.0j},
            vertex,
            reverse_wall,
            mode="no_prune",
        )
        self.assertTrue(reverse, "assembled node-H reverse column unexpectedly vanished")
        for output, amplitude in reverse.items():
            self.assertTrue(
                math.isfinite(amplitude.real) and math.isfinite(amplitude.imag),
                msg=f"non-finite reverse H amplitude at {output}",
            )
            if output in row:
                matrix[row[output], 1] = amplitude
        omitted = [
            amplitude for output, amplitude in reverse.items() if output not in row
        ]
        leakage.append({
            "count": len(omitted),
            "l2": math.sqrt(math.fsum(abs(value) ** 2 for value in omitted)),
            "max_abs": max((abs(value) for value in omitted), default=0.0),
        })

        self.assertEqual(matrix.shape, (2, 2))
        self.assertGreater(
            abs(matrix[1, 0]),
            0.0,
            "selected projected transition must be explicitly nonzero",
        )
        self.assertGreater(
            float(np.linalg.norm(matrix)),
            0.0,
            "projected matrix must not be identically zero",
        )
        defect = float(np.linalg.norm(matrix - matrix.conj().T))
        scale = max(float(np.linalg.norm(matrix)), np.finfo(float).tiny)
        relative_defect = defect / scale
        np.testing.assert_allclose(
            matrix,
            matrix.conj().T,
            rtol=2e-10,
            atol=2e-11,
            err_msg="non-vacuous projected production node-H matrix is not Hermitian",
        )
        print(
            "ASSEMBLED_NODE_H_NONVACUOUS_PROJECTED_BLOCK",
            {
                "vertex": vertex,
                "basis_dimension": len(basis),
                "input_max_doubled_spin": max(input_key[0]),
                "output_max_doubled_spin": max(output_key[0]),
                "forward_spin_wall_Jmax2": forward_wall,
                "reverse_spin_wall_Jmax2": reverse_wall,
                "selected_forward_amplitude": [
                    selected_amplitude.real, selected_amplitude.imag
                ],
                "matrix": [
                    [[z.real, z.imag] for z in line] for line in matrix
                ],
                "matrix_norm": float(np.linalg.norm(matrix)),
                "absolute_hermiticity_defect": defect,
                "relative_hermiticity_defect": relative_defect,
                "leakage_outside_projected_basis": leakage,
            },
        )
        self.assertEqual(len(leakage), 2)
        for item in leakage:
            self.assertGreaterEqual(item["count"], 0)
            self.assertGreaterEqual(item["l2"], 0.0)
            self.assertGreaterEqual(item["max_abs"], 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
