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


if __name__ == "__main__":
    unittest.main()
