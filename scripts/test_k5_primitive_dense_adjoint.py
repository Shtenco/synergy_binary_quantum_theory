#!/usr/bin/env python3
"""Independent finite-sector matrix adjoint oracle for a production K5 primitive.

The oracle materializes the complete closed-loop three-hit primitive and its reversed
production primitive on explicit orthonormal recoupling bases.  The finite
spin wall is part of the test definition; this is not a claim about an
unrestricted operator on the full Hilbert space.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

import k5_peter_weyl_safe_hda_column as k5


def basis_for_spins(spins):
    local_channels = [
        k5.allowed_k2_t(*k5.local_spins(spins, v))
        for v in k5.VERT
    ]
    return [
        (spins, tuple(channels))
        for channels in itertools.product(*local_channels)
    ]


def materialize_primitive(domain_in, domain_out, sequence, spin_wall):
    row = {key: i for i, key in enumerate(domain_out)}
    matrix = np.zeros((len(domain_out), len(domain_in)), dtype=complex)
    for col, key in enumerate(domain_in):
        result = k5.apply_primitive_oriented(
            key, sequence, spin_wall, mode="no_prune"
        )
        for out_key, amplitude in result.items():
            if out_key not in row:
                raise AssertionError(
                    f"production primitive escaped declared output basis: {out_key}"
                )
            matrix[row[out_key], col] += amplitude
    return matrix


class K5PrimitiveDenseAdjointTests(unittest.TestCase):
    def test_closed_triangle_primitive_matches_independent_reverse_matrix(self):
        # All ten links start at doubled spin 1 (j=1/2). The three edges of
        # the closed triangle may transition to doubled spin 0 or 2; Jmax2=2
        # excludes spin 3 in the reverse action.
        spin_wall = 2
        spins_in = (1,) * len(k5.EDGES)
        domain_in = basis_for_spins(spins_in)
        triangle_edges = tuple(k5.EIDX[e] for e in ((0, 1), (1, 2), (0, 2)))

        for iout, iin in itertools.product(range(2), repeat=2):
            with self.subTest(iout=iout, iin=iin):
                domain_out = []
                for edge_spins in itertools.product((0, 2), repeat=3):
                    spins_out = list(spins_in)
                    for edge_index, edge_spin in zip(triangle_edges, edge_spins):
                        spins_out[edge_index] = edge_spin
                    domain_out.extend(basis_for_spins(tuple(spins_out)))

                forward_path = (0, 1, 2, 0)
                forward_sequence = (("P", forward_path, iout, iin),)
                reverse_sequence = k5.adjoint_sequence(forward_sequence)

                forward = materialize_primitive(
                    domain_in, domain_out, forward_sequence, spin_wall
                )
                reverse = materialize_primitive(
                    domain_out, domain_in, reverse_sequence, spin_wall
                )

                np.testing.assert_allclose(
                    reverse,
                    forward.conj().T,
                    rtol=3e-12,
                    atol=3e-12,
                    err_msg=(
                        "finite-sector production primitive and independently "
                        "materialized reversed primitive are not adjoints"
                    ),
                )


    def test_local_recoupling_bases_are_orthonormal_on_all_test_sectors(self):
        # The dense adjoint comparison is meaningful only in orthonormal
        # domain/codomain bases. Check every local channel basis used by the
        # input sector and all eight output spin assignments explicitly.
        spin_configs = {(1,) * len(k5.EDGES)}
        triangle_edges = tuple(
            k5.EIDX[e] for e in ((0, 1), (1, 2), (0, 2))
        )
        for edge_spins in itertools.product((0, 2), repeat=3):
            spins = [1] * len(k5.EDGES)
            for edge_index, edge_spin in zip(triangle_edges, edge_spins):
                spins[edge_index] = edge_spin
            spin_configs.add(tuple(spins))

        for spins in sorted(spin_configs):
            for vertex in k5.VERT:
                local = k5.local_spins(spins, vertex)
                channels = k5.allowed_k2_t(*local)
                tensors = [
                    k5.oriented_intertwiner(vertex, local, channel)
                    for channel in channels
                ]
                gram = np.array([
                    [np.vdot(left, right) for right in tensors]
                    for left in tensors
                ], dtype=complex)
                np.testing.assert_allclose(
                    gram,
                    np.eye(len(channels)),
                    rtol=2e-13,
                    atol=2e-13,
                    err_msg=(
                        f"non-orthonormal local recoupling basis at "
                        f"vertex={vertex}, spins={local}"
                    ),
                )


    def test_all_twelve_oriented_specs_have_exact_combinatorics(self):
        # Four choices of omitted neighbor, each with three cyclic orders.
        # The sign is fixed by the omitted-neighbor position, not by the
        # cyclic rotation. Compare the generator to an independently written
        # expected list and also enforce uniqueness/completeness.
        expected = []
        for omitted in range(4):
            triple = tuple(w for w in k5.NEIG[0] if w != k5.NEIG[0][omitted])
            a, b, c = triple
            sign = (-1) ** omitted
            expected.extend([
                (sign, (0, a, b, c)),
                (sign, (0, b, c, a)),
                (sign, (0, c, a, b)),
            ])

        actual = k5.oriented_specs(0)
        self.assertEqual(len(actual), 12)
        self.assertEqual(actual, tuple(expected))
        self.assertEqual(len(set(actual)), 12)

        # The same combinatorial rule must hold at every K5 vertex.
        for vertex in k5.VERT:
            specs = k5.oriented_specs(vertex)
            self.assertEqual(len(specs), 12)
            self.assertEqual(len(set(specs)), 12)
            neigh = k5.NEIG[vertex]
            independently_expected = []
            for omitted in range(4):
                triple = tuple(w for w in neigh if w != neigh[omitted])
                a, b, c = triple
                sign = (-1) ** omitted
                independently_expected.extend([
                    (sign, (vertex, a, b, c)),
                    (sign, (vertex, b, c, a)),
                    (sign, (vertex, c, a, b)),
                ])
            self.assertEqual(specs, tuple(independently_expected))

    def test_full_32_term_T_sum_matches_reverse_matrix_element(self):
        # Test the coherent T sum itself, not merely its 32 monomial pieces.
        # Accumulate the forward and reverse matrix elements independently
        # with fsum, preserving the coefficient signs and cross-term addition.
        import math
        from collections import defaultdict

        vertex, a, b, c = 0, 1, 2, 3
        spin_wall = 2
        input_key = ((1,) * len(k5.EDGES), (0,) * len(k5.VERT))
        terms = k5.T_sequences(vertex, a, b, c)
        self.assertEqual(len(terms), 32)

        forward_terms = defaultdict(list)
        for coefficient, sequence in terms:
            column = k5.apply_primitive_oriented(
                input_key, sequence, spin_wall, mode="no_prune"
            )
            for output_key, amplitude in column.items():
                forward_terms[output_key].append(coefficient * amplitude)

        forward = {
            key: complex(
                math.fsum(value.real for value in values),
                math.fsum(value.imag for value in values),
            )
            for key, values in forward_terms.items()
        }
        forward = {key: value for key, value in forward.items() if value != 0j}
        self.assertTrue(forward, "the complete 32-term T column unexpectedly vanished")

        # Select the largest-magnitude nonzero element deterministically.
        output_key = min(
            forward,
            key=lambda key: (-abs(forward[key]), repr(key)),
        )
        forward_element = forward[output_key]

        reverse_contributions = []
        for coefficient, sequence in terms:
            reverse_sequence = k5.adjoint_sequence(sequence)
            reverse_column = k5.apply_primitive_oriented(
                output_key, reverse_sequence, spin_wall, mode="no_prune"
            )
            reverse_contributions.append(
                coefficient * reverse_column.get(input_key, 0j)
            )
        reverse_element = complex(
            math.fsum(value.real for value in reverse_contributions),
            math.fsum(value.imag for value in reverse_contributions),
        )

        np.testing.assert_allclose(
            reverse_element,
            np.conjugate(forward_element),
            rtol=1e-10,
            atol=1e-11,
            err_msg=(
                "the coherent 32-term production T sum disagrees with "
                "the independently accumulated reversed matrix element"
            ),
        )

    def test_full_five_hit_T_monomial_matches_reverse_matrix_elements(self):
        # One complete production monomial from T_sequences contains the
        # volume insertion and all three holonomy paths (five link hits).
        # Rather than infer its adjoint from sequence reversal alone, compare
        # every nonzero forward matrix element from a fixed basis ket with
        # the independently evaluated reversed production action.
        vertex, a, b, c = 0, 1, 2, 3
        spin_wall = 2
        input_key = ((1,) * len(k5.EDGES), (0,) * len(k5.VERT))
        terms = k5.T_sequences(vertex, a, b, c)
        self.assertEqual(len(terms), 32)
        coefficient, sequence = terms[0]
        self.assertEqual(coefficient, 1)

        forward = k5.apply_primitive_oriented(
            input_key, sequence, spin_wall, mode="no_prune"
        )
        self.assertTrue(forward, "selected complete T monomial unexpectedly vanished")
        reverse_sequence = k5.adjoint_sequence(sequence)

        for output_key, amplitude in forward.items():
            reverse = k5.apply_primitive_oriented(
                output_key, reverse_sequence, spin_wall, mode="no_prune"
            )
            reverse_matrix_element = reverse.get(input_key, 0j)
            with self.subTest(output_key=output_key):
                np.testing.assert_allclose(
                    reverse_matrix_element,
                    np.conjugate(amplitude),
                    rtol=8e-11,
                    atol=8e-12,
                    err_msg=(
                        "complete five-hit T monomial disagrees with its "
                        "reversed production matrix element"
                    ),
                )


    def test_full_32_term_T_two_by_two_matrix_block_matches_reverse(self):
        # A stronger finite check than a single selected matrix element:
        # build two orthonormal input columns, choose two output rows from
        # their coherent forward support, and independently build the reverse
        # T^dagger columns for those rows. This is a projected 2x2 block, not
        # a claim that either selected subspace is invariant.
        import math
        from collections import defaultdict

        vertex, a, b, c = 0, 1, 2, 3
        spin_wall = 2
        spins = (1,) * len(k5.EDGES)
        input_keys = (
            (spins, (0,) * len(k5.VERT)),
            (spins, (0, 0, 0, 0, 2)),
        )
        for key in input_keys:
            for v, channel in enumerate(key[1]):
                self.assertIn(
                    channel, k5.allowed_k2_t(*k5.local_spins(key[0], v)),
                    msg=f"input channel {channel} is not allowed at vertex {v}",
                )

        terms = k5.T_sequences(vertex, a, b, c)
        self.assertEqual(len(terms), 32)

        def coherent_column(input_key, reverse=False):
            contributions = defaultdict(list)
            for coefficient, sequence in terms:
                action = k5.adjoint_sequence(sequence) if reverse else sequence
                mapping = k5.apply_primitive_oriented(
                    input_key, action, spin_wall, mode="no_prune"
                )
                for output_key, amplitude in mapping.items():
                    contributions[output_key].append(coefficient * amplitude)
            result = {}
            for output_key, values in contributions.items():
                amplitude = complex(
                    math.fsum(value.real for value in values),
                    math.fsum(value.imag for value in values),
                )
                if amplitude != 0j:
                    result[output_key] = amplitude
            return result

        forward_columns = [coherent_column(key) for key in input_keys]
        union_support = set().union(*(set(column) for column in forward_columns))
        self.assertGreaterEqual(
            len(union_support), 2,
            "the two input columns do not produce two distinct output rows",
        )

        # Deterministic row choice: largest maximum column amplitude, then
        # repr(key) as a stable tie-breaker.
        output_keys = tuple(sorted(
            union_support,
            key=lambda key: (
                -max(abs(column.get(key, 0j)) for column in forward_columns),
                repr(key),
            ),
        )[:2])

        forward_block = np.array([
            [column.get(output_key, 0j) for column in forward_columns]
            for output_key in output_keys
        ], dtype=complex)

        reverse_columns = [coherent_column(key, reverse=True) for key in output_keys]
        reverse_block = np.array([
            [column.get(input_key, 0j) for column in reverse_columns]
            for input_key in input_keys
        ], dtype=complex)

        # Explicitly retain the projected-block boundary ledger. These are
        # amplitudes outside the selected 2x2 block, not transitions clipped
        # by the spin wall; the latter require a separate unbounded comparison.
        forward_leakage = []
        for column in forward_columns:
            omitted = [
                amplitude for key, amplitude in column.items()
                if key not in output_keys
            ]
            forward_leakage.append({
                "count": len(omitted),
                "l2": math.sqrt(math.fsum(abs(value) ** 2 for value in omitted)),
                "max_abs": max((abs(value) for value in omitted), default=0.0),
            })

        reverse_leakage = []
        for column in reverse_columns:
            omitted = [
                amplitude for key, amplitude in column.items()
                if key not in input_keys
            ]
            reverse_leakage.append({
                "count": len(omitted),
                "l2": math.sqrt(math.fsum(abs(value) ** 2 for value in omitted)),
                "max_abs": max((abs(value) for value in omitted), default=0.0),
            })

        self.assertEqual(forward_block.shape, (2, 2))
        self.assertEqual(reverse_block.shape, (2, 2))
        np.testing.assert_allclose(
            reverse_block,
            forward_block.conj().T,
            rtol=1e-10,
            atol=1e-11,
            err_msg=(
                "the independently accumulated reverse 32-term T block "
                "does not equal the conjugate transpose of the forward block"
            ),
        )

        # The ledger is part of the test's output contract: omitted support
        # is measured and never silently interpreted as a closed subspace.
        self.assertEqual(len(forward_leakage), 2)
        self.assertEqual(len(reverse_leakage), 2)
        print(
            "T_BLOCK_LEAKAGE",
            {
                "forward": forward_leakage,
                "reverse": reverse_leakage,
                "spin_wall_Jmax2": spin_wall,
                "block_shape": (2, 2),
            },
        )
        for item in forward_leakage + reverse_leakage:
            self.assertGreaterEqual(item["count"], 0)
            self.assertGreaterEqual(item["l2"], 0.0)
            self.assertGreaterEqual(item["max_abs"], 0.0)


    def test_full_32_term_T_spin_wall_2_vs_4_sensitivity_ledger(self):
        # The primitive hit rule changes an edge spin by +/-1 (doubled-spin
        # units) at every hit. Compare two finite walls on the same coherent
        # T column to expose paths that are clipped at Jmax2=2 but survive at
        # Jmax2=4. This is a wall-sensitivity diagnostic, not an unbounded
        # reference and not a proof that Jmax2=4 is sufficient.
        import math
        from collections import defaultdict

        vertex, a, b, c = 0, 1, 2, 3
        input_key = ((1,) * len(k5.EDGES), (0,) * len(k5.VERT))
        terms = k5.T_sequences(vertex, a, b, c)
        self.assertEqual(len(terms), 32)

        def coherent_at_wall(wall):
            contributions = defaultdict(list)
            for coefficient, sequence in terms:
                mapping = k5.apply_primitive_oriented(
                    input_key, sequence, wall, mode="no_prune"
                )
                for output_key, amplitude in mapping.items():
                    contributions[output_key].append(coefficient * amplitude)
            result = {}
            for output_key, values in contributions.items():
                amplitude = complex(
                    math.fsum(value.real for value in values),
                    math.fsum(value.imag for value in values),
                )
                if amplitude != 0j:
                    result[output_key] = amplitude
            return result

        wall2 = coherent_at_wall(2)
        wall4 = coherent_at_wall(4)
        wall6 = coherent_at_wall(6)
        self.assertTrue(wall2, "T column vanished at Jmax2=2")
        self.assertTrue(wall4, "T column vanished at Jmax2=4")
        self.assertTrue(wall6, "T column vanished at Jmax2=6")
        # For this exact input and each T monomial, the c-v edge is hit at
        # most twice and every other edge at most once. Starting from doubled
        # spin 1, no reachable edge can exceed doubled spin 3. Therefore both
        # walls 4 and 6 are above the analytically derived reachable ceiling;
        # identical maps certify wall saturation for this one input sector.
        self.assertEqual(wall4, wall6)
        for label, mapping in (("wall2", wall2), ("wall4", wall4), ("wall6", wall6)):
            for key, value in mapping.items():
                self.assertTrue(
                    np.isfinite(value.real) and np.isfinite(value.imag),
                    msg=f"non-finite amplitude in {label} at {key}",
                )

        common = set(wall2) & set(wall4)
        only_wall2 = set(wall2) - set(wall4)
        only_wall4 = set(wall4) - set(wall2)
        low_wall_delta = [
            wall4[key] - wall2[key]
            for key in common
            if max(key[0]) <= 2
        ]
        low_l2 = math.sqrt(math.fsum(abs(value) ** 2 for value in low_wall_delta))
        low_max = max((abs(value) for value in low_wall_delta), default=0.0)
        high_wall_keys = [
            key for key in wall4 if max(key[0]) > 2
        ]
        high_wall_l2 = math.sqrt(math.fsum(
            abs(wall4[key]) ** 2 for key in high_wall_keys
        ))
        print(
            "T_SPIN_WALL_SENSITIVITY",
            {
                "input": input_key,
                "walls": (2, 4, 6),
                "derived_max_reachable_doubled_spin": 3,
                "wall4_equals_wall6": wall4 == wall6,
                "support_wall2": len(wall2),
                "support_wall4": len(wall4),
                "common_support": len(common),
                "only_wall2": len(only_wall2),
                "only_wall4": len(only_wall4),
                "wall4_outputs_above_wall2": len(high_wall_keys),
                "wall4_high_spin_l2": high_wall_l2,
                "common_low_spin_delta_l2": low_l2,
                "common_low_spin_delta_max": low_max,
            },
        )
        self.assertGreaterEqual(low_l2, 0.0)
        self.assertGreaterEqual(high_wall_l2, 0.0)



    def test_all_k5_T_and_H_monomials_have_edge_hit_multiplicity_at_most_two(self):
        # A path hit changes one edge's doubled spin by exactly +/-1.
        # Every T monomial has two one-hit paths on (c,v) plus one
        # three-hit closed path (v,a,b,v). Verify this combinatorially for
        # every vertex, all 12 oriented specs, and all 32 signed monomials.
        # Reversal preserves undirected edge multiplicities, so the same
        # bound applies to T^dagger and to each summand of H_v.
        from collections import Counter

        checked = 0
        multiplicity_histogram = Counter()
        for vertex in k5.VERT:
            specs = k5.oriented_specs(vertex)
            self.assertEqual(len(specs), 12)
            for sign, (v, a, b, c) in specs:
                self.assertEqual(v, vertex)
                self.assertEqual(len({v, a, b, c}), 4)
                terms = k5.T_sequences(v, a, b, c)
                self.assertEqual(len(terms), 32)
                for coefficient, sequence in terms:
                    self.assertIn(coefficient, (-1, 1))
                    edge_hits = Counter()
                    total_hits = 0
                    for operation in sequence:
                        if operation[0] != "P":
                            self.assertEqual(operation[0], "V")
                            continue
                        path = operation[1]
                        for x, y in zip(path[:-1], path[1:]):
                            edge_hits[tuple(sorted((x, y)))] += 1
                            total_hits += 1
                    self.assertEqual(total_hits, 5)
                    self.assertEqual(len(edge_hits), 4)
                    self.assertEqual(sorted(edge_hits.values()), [1, 1, 1, 2])
                    self.assertLessEqual(max(edge_hits.values()), 2)
                    multiplicity_histogram[tuple(sorted(edge_hits.values()))] += 1
                    # Reversing a path/sequence preserves the undirected
                    # edge-hit multiset, so this same bound covers T^dagger.
                    reversed_sequence = k5.adjoint_sequence(sequence)
                    reverse_hits = Counter()
                    for operation in reversed_sequence:
                        if operation[0] == "P":
                            for x, y in zip(operation[1][:-1], operation[1][1:]):
                                reverse_hits[tuple(sorted((x, y)))] += 1
                    self.assertEqual(reverse_hits, edge_hits)
                    checked += 1

        self.assertEqual(checked, 5 * 12 * 32)
        self.assertEqual(multiplicity_histogram[(1, 1, 1, 2)], checked)


if __name__ == "__main__":
    unittest.main()
