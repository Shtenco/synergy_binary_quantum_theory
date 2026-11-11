#!/usr/bin/env python3
"""Non-vacuous finite-sector Hermiticity oracle for the assembled K5 node H.

The second basis ket is selected deterministically from a nonzero production
H action on the first ket. The complete projected 2x2 matrix contains a
measured nonzero transition. The reverse matrix element is computed
independently from only those signed production specifications whose exact
link-hit support can connect the selected pair; all omitted specifications
are proven unable to contribute to that matrix element. The diagonal entries
vanish by the independently tested (1,1,1,2) hit-multiplicity selection rule.
This is a two-state finite-sector test, not a full-sector or continuum theorem.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

import k5_peter_weyl_safe_hda_column as k5


def sequence_hit_edges(sequence):
    return {
        tuple(sorted((x, y)))
        for op in sequence
        if op[0] == "P"
        for x, y in zip(op[1][:-1], op[1][1:])
    }


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

        # Materialize the full no-prune H column once. S_in=1, so the
        # universal per-edge hit bound gives the exact sufficient wall 3.
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
        output_key, selected_amplitude = sorted(
            candidates, key=lambda item: (-abs(item[1]), item[0])
        )[0]
        validate_key(output_key)

        basis = (input_key, output_key)
        row = {key: i for i, key in enumerate(basis)}
        matrix = np.zeros((2, 2), dtype=complex)
        matrix[0, 0] = forward.get(input_key, 0j)
        matrix[1, 0] = forward[output_key]

        forward_omitted = [
            amplitude for output, amplitude in forward.items() if output not in row
        ]
        forward_leakage = {
            "count": len(forward_omitted),
            "l2": math.sqrt(math.fsum(abs(value) ** 2 for value in forward_omitted)),
            "max_abs": max((abs(value) for value in forward_omitted), default=0.0),
        }

        # For a fixed pair, a specification can contribute to <input|H|output>
        # only if every edge whose spin differs is in that specification's
        # hit support. This is an exact support filter, not amplitude pruning.
        changed_edges = {
            edge
            for edge, sin, sout in zip(k5.EDGES, output_key[0], input_key[0])
            if sin != sout
        }
        self.assertIn(
            len(changed_edges),
            (3, 4),
            "a five-hit monomial must change its three single-hit edges and may change the double-hit edge",
        )
        reverse_terms = []
        candidate_specs = []
        omitted_specs = []
        reverse_wall = max(output_key[0]) + 2
        for sign, spec in k5.oriented_specs(vertex):
            sequences = k5.T_sequences(*spec)
            hit_support = sequence_hit_edges(sequences[0][1])
            self.assertTrue(
                all(sequence_hit_edges(seq) == hit_support for _, seq in sequences),
                msg=f"monomial hit support unexpectedly varies within spec {spec}",
            )
            if not changed_edges.issubset(hit_support):
                omitted_specs.append(spec)
                continue

            candidate_specs.append(spec)
            t_map = k5.apply_T_cached_state(
                {output_key: 1.0 + 0.0j},
                spec,
                reverse_wall,
                adj=False,
                mode="no_prune",
            )
            td_map = k5.apply_T_cached_state(
                {output_key: 1.0 + 0.0j},
                spec,
                reverse_wall,
                adj=True,
                mode="no_prune",
            )
            reverse_terms.extend((
                0.5 * sign * t_map.get(input_key, 0j),
                0.5 * sign * td_map.get(input_key, 0j),
            ))

        reverse_amplitude = complex(
            math.fsum(value.real for value in reverse_terms),
            math.fsum(value.imag for value in reverse_terms),
        )
        matrix[0, 1] = reverse_amplitude

        # Every monomial has three single-hit edges, so its link-spin label
        # vector cannot be unchanged. Both projected diagonal entries are
        # exactly zero by the committed exhaustive combinatorial test.
        self.assertEqual(matrix[0, 0], 0j)
        matrix[1, 1] = 0j

        self.assertGreater(abs(matrix[1, 0]), 0.0)
        self.assertGreater(float(np.linalg.norm(matrix)), 0.0)
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
                "reverse_amplitude": [
                    reverse_amplitude.real, reverse_amplitude.imag
                ],
                "matrix": [
                    [[z.real, z.imag] for z in line] for line in matrix
                ],
                "matrix_norm": float(np.linalg.norm(matrix)),
                "absolute_hermiticity_defect": defect,
                "relative_hermiticity_defect": relative_defect,
                "forward_leakage_outside_projected_basis": forward_leakage,
                "reverse_candidate_spec_count": len(candidate_specs),
                "reverse_omitted_spec_count_proven_zero_for_pair": len(omitted_specs),
            },
        )
        self.assertGreater(len(candidate_specs), 0)
        self.assertEqual(
            len(candidate_specs) + len(omitted_specs),
            len(k5.oriented_specs(vertex)),
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
