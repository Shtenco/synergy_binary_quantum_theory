#!/usr/bin/env python3
from __future__ import annotations

import functools

import numpy as np

import bqg_safe_boundary_master_gate as R
import bqg_covariant_volume_recoupling as Rec


@functools.lru_cache(None)
def oriented_coupled_tensor(v, spins_local, K, L, J, M):
    tensor = Rec.coupled_tensor(tuple(spins_local), K, L, J, M).copy()
    for leg, neighbor in enumerate(R.NEIG[v]):
        if neighbor < v:
            tensor = R.apply_axis_np(
                tensor, leg, R.epsilon_j(spins_local[leg])
            )
    return tensor


@functools.lru_cache(None)
def local_hit_outputs(
    v, spins_local, leg, s_out, kind, color, label, tol=1e-12
):
    """Reduced one-leg transfer; tol=0 disables magnitude pruning."""
    spins_local = tuple(spins_local)
    K, L, J, M = label
    s_in = spins_local[leg]
    magnetic_component = 1 if color == 0 else -1

    if kind == "cg":
        matrix = R.cg_leg_matrix(s_in, s_out, magnetic_component)
    elif kind == "adj":
        matrix = R.cg_leg_matrix(s_out, s_in, magnetic_component).conj().T
    else:
        raise ValueError(kind)

    output_spins = list(spins_local)
    output_spins[leg] = s_out
    output_spins = tuple(output_spins)
    basis_in = oriented_coupled_tensor(v, spins_local, K, L, J, M)
    transformed = R.apply_axis_np(basis_in, leg, matrix)

    outputs = []
    for J_out in Rec.allowed_total_J(output_spins):
        if abs(J_out - J) != 1:
            continue
        for M_out in R.m2vals_t(J_out):
            for K_out, L_out in Rec.multiplicity_labels(output_spins, J_out):
                basis_out = oriented_coupled_tensor(
                    v, output_spins, K_out, L_out, J_out, M_out
                )
                amplitude = np.vdot(basis_out, transformed)
                if abs(amplitude) > tol:
                    outputs.append(
                        ((K_out, L_out, J_out, M_out), complex(amplitude))
                    )
    return tuple(outputs)
