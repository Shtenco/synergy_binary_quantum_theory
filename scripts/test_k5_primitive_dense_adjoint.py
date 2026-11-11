#!/usr/bin/env python3
"""Independent finite-sector matrix adjoint oracle for a production K5 primitive.

The oracle materializes the complete one-link-hit primitive and its reversed
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
    def test_one_link_hit_primitive_matches_independent_reverse_matrix(self):
        # All ten links start at doubled spin 1 (j=1/2).  Only edge (0,1)
        # may transition to doubled spin 0 or 2; Jmax2=2 excludes spin 3
        # in the reverse action, making this an explicitly finite-wall test.
        spin_wall = 2
        spins_in = (1,) * len(k5.EDGES)
        domain_in = basis_for_spins(spins_in)
        edge_index = k5.EIDX[(0, 1)]

        for iout, iin in itertools.product(range(2), repeat=2):
            with self.subTest(iout=iout, iin=iin):
                domain_out = []
                for edge_spin in (0, 2):
                    spins_out = list(spins_in)
                    spins_out[edge_index] = edge_spin
                    domain_out.extend(basis_for_spins(tuple(spins_out)))

                forward_path = (0, 1)
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


if __name__ == "__main__":
    unittest.main()
