#!/usr/bin/env python3
"""Rectangular projected finite-sector adjointness gate for the assembled K5 node H.

Two orthonormal K-channel inputs share the all-link-spin-1 carrier. Their
production no-prune H columns define a union of outputs; a deterministic top-k
subset forms a rectangular projected block. Reverse entries are recomputed
independently from production T/T† maps after an exact link-hit-support filter.
Omitted forward support is measured per input column. This is not a complete
finite-sector or continuum self-adjointness theorem.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

import k5_peter_weyl_safe_hda_column as k5


TOP_K = 4


def sequence_hit_edges(sequence):
    return {
        tuple(sorted((x, y)))
        for op in sequence
        if op[0] == "P"
        for x, y in zip(op[1][:-1], op[1][1:])
    }


def stable_complex_sum(values):
    return complex(
        math.fsum(value.real for value in values),
        math.fsum(value.imag for value in values),
    )


class K5RectangularNodeHAdjointnessTests(unittest.TestCase):
    def test_rectangular_projected_node_H_block_is_adjoint(self):
        vertex = 0
        spins = (1,) * len(k5.EDGES)
        inputs = (
            (spins, (0, 0, 0, 0, 0)),
            (spins, (2, 0, 0, 0, 0)),
        )
        self.assertEqual(len(set(inputs)), len(inputs))

        for key in inputs:
            for v, channel in enumerate(key[1]):
                self.assertIn(
                    channel,
                    k5.allowed_k2_t(*k5.local_spins(key[0], v)),
                    msg=f"invalid recoupling channel {channel} at vertex {v}",
                )

        forward_wall = max(max(key[0]) for key in inputs) + 2
        forward_columns = []
        for key in inputs:
            column = k5.apply_H_cached_state(
                {key: 1.0 + 0.0j},
                vertex,
                forward_wall,
                mode="no_prune",
            )
            self.assertTrue(column, f"assembled H column vanished for {key}")
            for output_key, amplitude in column.items():
                self.assertTrue(
                    math.isfinite(amplitude.real) and math.isfinite(amplitude.imag),
                    msg=f"non-finite forward amplitude at {output_key}",
                )
            forward_columns.append(column)

        # Rank the union by its Euclidean column-amplitude score; lexicographic
        # state ordering makes ties reproducible across runs.
        reachable = set().union(*(column.keys() for column in forward_columns))
        self.assertGreaterEqual(len(reachable), TOP_K)
        ranked = sorted(
            reachable,
            key=lambda key: (
                -math.sqrt(math.fsum(abs(column.get(key, 0j)) ** 2
                                     for column in forward_columns)),
                key,
            ),
        )
        # Seed one nonzero selected output per input column before filling
        # the remaining slots globally. This prevents a nominal 4x2 block
        # from silently containing an all-zero input column.
        outputs_list = []
        for column in forward_columns:
            best_for_column = min(
                (key for key in reachable if column.get(key, 0j) != 0j),
                key=lambda key: (-abs(column[key]), key),
            )
            if best_for_column not in outputs_list:
                outputs_list.append(best_for_column)
        for key in ranked:
            if key not in outputs_list:
                outputs_list.append(key)
            if len(outputs_list) == TOP_K:
                break
        outputs = tuple(outputs_list)
        self.assertEqual(len(set(outputs)), TOP_K)
        self.assertTrue(
            all(any(column.get(key, 0j) != 0j for key in outputs)
                for column in forward_columns),
            "each input column must contribute a nonzero selected matrix element",
        )

        # A[y,x] is taken directly from complete no-prune forward columns.
        A = np.zeros((len(outputs), len(inputs)), dtype=complex)
        for iy, output_key in enumerate(outputs):
            for ix, column in enumerate(forward_columns):
                A[iy, ix] = column.get(output_key, 0j)

        leakage = []
        selected = set(outputs)
        for column in forward_columns:
            omitted = [amp for key, amp in column.items() if key not in selected]
            leakage.append({
                "count": len(omitted),
                "l1": math.fsum(abs(value) for value in omitted),
                "l2": math.sqrt(math.fsum(abs(value) ** 2 for value in omitted)),
                "max_abs": max((abs(value) for value in omitted), default=0.0),
            })

        # Both inputs have the same link-spin vector. Thus the exact changed
        # edge set and reverse support filter can be shared across columns.
        B = np.zeros((len(inputs), len(outputs)), dtype=complex)
        candidate_counts = []
        omitted_spec_counts = []
        for iy, output_key in enumerate(outputs):
            changed_edges = {
                edge
                for edge, sin, sout in zip(k5.EDGES, output_key[0], spins)
                if sin != sout
            }
            self.assertIn(
                len(changed_edges),
                (3, 4),
                "every five-hit monomial changes its three single-hit edges",
            )

            reverse_wall = max(output_key[0]) + 2
            reverse_terms = [[] for _ in inputs]
            candidate_specs = []
            omitted_specs = []
            for sign, spec in k5.oriented_specs(vertex):
                sequences = k5.T_sequences(*spec)
                hit_support = sequence_hit_edges(sequences[0][1])
                self.assertTrue(
                    all(sequence_hit_edges(seq) == hit_support
                        for _, seq in sequences),
                    msg=f"monomial support varies within oriented spec {spec}",
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
                for ix, input_key in enumerate(inputs):
                    reverse_terms[ix].extend((
                        0.5 * sign * t_map.get(input_key, 0j),
                        0.5 * sign * td_map.get(input_key, 0j),
                    ))

            for ix in range(len(inputs)):
                B[ix, iy] = stable_complex_sum(reverse_terms[ix])
            candidate_counts.append(len(candidate_specs))
            omitted_spec_counts.append(len(omitted_specs))

        self.assertTrue(np.any(np.abs(A) > 0.0), "projected forward block is vacuous")
        self.assertTrue(np.all(np.isfinite(A.real)) and np.all(np.isfinite(A.imag)))
        self.assertTrue(np.all(np.isfinite(B.real)) and np.all(np.isfinite(B.imag)))

        residual = A - B.conj().T
        frobenius_defect = float(np.linalg.norm(residual))
        max_entry_defect = float(np.max(np.abs(residual)))
        relative_frobenius_defect = frobenius_defect / max(
            float(np.linalg.norm(A)), float(np.linalg.norm(B)), np.finfo(float).tiny
        )
        np.testing.assert_allclose(
            A,
            B.conj().T,
            rtol=2e-10,
            atol=2e-11,
            err_msg="rectangular projected production node-H block is not adjoint",
        )

        print(
            "RECTANGULAR_PROJECTED_NODE_H_ADJOINTNESS",
            {
                "vertex": vertex,
                "input_dimension": len(inputs),
                "selected_output_dimension": len(outputs),
                "complete_reachable_output_count": len(reachable),
                "top_k": TOP_K,
                "forward_spin_wall_Jmax2": forward_wall,
                "reverse_spin_walls_Jmax2": [
                    max(output_key[0]) + 2 for output_key in outputs
                ],
                "forward_leakage_per_input_column": leakage,
                "reverse_candidate_specs_per_output": candidate_counts,
                "reverse_omitted_specs_proven_zero_per_output": omitted_spec_counts,
                "forward_matrix_norm": float(np.linalg.norm(A)),
                "reverse_matrix_norm": float(np.linalg.norm(B)),
                "frobenius_adjoint_defect": frobenius_defect,
                "max_entry_adjoint_defect": max_entry_defect,
                "relative_frobenius_adjoint_defect": relative_frobenius_defect,
                "row_residual_norms": [
                    float(np.linalg.norm(residual[i, :]))
                    for i in range(residual.shape[0])
                ],
                "column_residual_norms": [
                    float(np.linalg.norm(residual[:, j]))
                    for j in range(residual.shape[1])
                ],
                "scope": "top-k rectangular projection; omitted forward leakage reported",
            },
        )
        self.assertEqual(len(candidate_counts), len(outputs))
        self.assertTrue(all(0 < n <= len(k5.oriented_specs(vertex))
                            for n in candidate_counts))
        self.assertEqual(
            [c + o for c, o in zip(candidate_counts, omitted_spec_counts)],
            [len(k5.oriented_specs(vertex))] * len(outputs),
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
