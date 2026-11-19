#!/usr/bin/env python3
"""Finite graph Dirichlet-form primitives for continuum-limit diagnostics.

This module computes exact finite-level energies, mass norms, generators,
resolvents, and prolongation/recovery probe diagnostics. Numerical probes can
falsify a proposed convergence scheme; they do not prove Mosco convergence.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GraphDirichletForm:
    n_vertices: int
    edges: np.ndarray
    edge_weights: np.ndarray
    vertex_measure: np.ndarray
    scale: float
    stiffness: object

    @classmethod
    def build(cls, n_vertices, edges, edge_weights, vertex_measure, *, scale=1.0):
        from scipy.sparse import coo_matrix

        n = int(n_vertices)
        e = np.asarray(edges)
        w = np.asarray(edge_weights, dtype=float)
        mu = np.asarray(vertex_measure, dtype=float)
        c = float(scale)
        if n < 1:
            raise ValueError("n_vertices must be positive")
        if e.ndim != 2 or e.shape[1] != 2 or not np.issubdtype(e.dtype, np.integer):
            raise ValueError("edges must be an integer array with shape (m, 2)")
        if w.shape != (len(e),) or not np.all(np.isfinite(w)) or np.any(w <= 0):
            raise ValueError("edge_weights must be finite, positive, and match edges")
        if mu.shape != (n,) or not np.all(np.isfinite(mu)) or np.any(mu <= 0):
            raise ValueError("vertex_measure must be finite, positive, and match vertices")
        if not np.isfinite(c) or c <= 0:
            raise ValueError("scale must be finite and positive")
        if len(e):
            if np.any(e < 0) or np.any(e >= n):
                raise ValueError("edge endpoint outside graph")
            if np.any(e[:, 0] == e[:, 1]):
                raise ValueError("self-loops are not supported")
            canonical = np.sort(e, axis=1)
            if len(np.unique(canonical, axis=0)) != len(e):
                raise ValueError("duplicate undirected edges are not supported")
        e = e.astype(np.int64, copy=False)
        rows = np.concatenate((e[:, 0], e[:, 1], e[:, 0], e[:, 1]))
        cols = np.concatenate((e[:, 0], e[:, 1], e[:, 1], e[:, 0]))
        data = np.concatenate((c * w, c * w, -c * w, -c * w))
        stiffness = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
        stiffness.sum_duplicates()
        return cls(n, e, w, mu, c, stiffness)

    def energy(self, values):
        u = np.asarray(values, dtype=float)
        if u.shape != (self.n_vertices,) or not np.all(np.isfinite(u)):
            raise ValueError("values must be a finite vector on the graph")
        differences = u[self.edges[:, 0]] - u[self.edges[:, 1]]
        return float(self.scale * np.dot(self.edge_weights, differences * differences))

    def mass_norm_squared(self, values):
        u = np.asarray(values, dtype=float)
        if u.shape != (self.n_vertices,) or not np.all(np.isfinite(u)):
            raise ValueError("values must be a finite vector on the graph")
        return float(np.dot(self.vertex_measure, u * u))

    def apply_generator(self, values):
        u = np.asarray(values, dtype=float)
        if u.shape != (self.n_vertices,) or not np.all(np.isfinite(u)):
            raise ValueError("values must be a finite vector on the graph")
        return np.asarray(self.stiffness @ u).ravel() / self.vertex_measure

    def resolvent_action(self, rhs, *, z):
        from scipy.sparse import diags
        from scipy.sparse.linalg import spsolve

        f = np.asarray(rhs, dtype=float)
        spectral_parameter = float(z)
        if f.shape != (self.n_vertices,) or not np.all(np.isfinite(f)):
            raise ValueError("rhs must be a finite vector on the graph")
        if not np.isfinite(spectral_parameter) or spectral_parameter <= 0:
            raise ValueError("z must be finite and positive")
        mass = diags(self.vertex_measure, format="csr")
        system = self.stiffness + spectral_parameter * mass
        solution = np.asarray(spsolve(system, self.vertex_measure * f)).ravel()
        residual = np.asarray(system @ solution).ravel() - self.vertex_measure * f
        relative_residual = float(
            np.linalg.norm(residual) /
            max(np.linalg.norm(self.vertex_measure * f), np.finfo(float).tiny)
        )
        if not np.all(np.isfinite(solution)) or not np.isfinite(relative_residual):
            raise RuntimeError("resolvent solve produced non-finite values")
        return solution, relative_residual

    def summary(self):
        return {
            "vertices": self.n_vertices,
            "edges": int(len(self.edges)),
            "measure_total": float(np.sum(self.vertex_measure)),
            "scale": self.scale,
            "stiffness_symmetry_error": float(
                np.max(np.abs((self.stiffness - self.stiffness.T).data), initial=0.0)
            ),
            "stiffness_row_sum_max_abs": float(
                np.max(np.abs(np.asarray(self.stiffness.sum(axis=1)).ravel()), initial=0.0)
            ),
            "form_convention": "E(u) = scale * sum_edges(weight * (u_i-u_j)^2)",
        }


def validate_prolongation(prolongation, *, n_fine, n_coarse, tolerance=1e-12):
    """Validate a nonnegative constant-preserving coarse-to-fine map."""
    import scipy.sparse as sp

    p = sp.csr_matrix(prolongation, dtype=float)
    if p.shape != (int(n_fine), int(n_coarse)):
        raise ValueError("prolongation has the wrong shape")
    if not np.all(np.isfinite(p.data)) or np.any(p.data < -tolerance):
        raise ValueError("prolongation must be finite and nonnegative")
    row_sums = np.asarray(p.sum(axis=1)).ravel()
    error = float(np.max(np.abs(row_sums - 1.0), initial=0.0))
    if error > tolerance:
        raise ValueError(f"prolongation does not preserve constants: {error}")
    return p, error


def recovery_probe_report(coarse_form, fine_form, prolongation, coarse_vectors):
    """Measure candidate recovery-sequence energies; this is not a Mosco proof."""
    p, constant_error = validate_prolongation(
        prolongation, n_fine=fine_form.n_vertices,
        n_coarse=coarse_form.n_vertices,
    )
    probes = []
    for index, vector in enumerate(coarse_vectors):
        coarse = np.asarray(vector, dtype=float)
        if coarse.shape != (coarse_form.n_vertices,):
            raise ValueError("coarse probe has the wrong shape")
        fine = np.asarray(p @ coarse).ravel()
        coarse_energy = coarse_form.energy(coarse)
        fine_energy = fine_form.energy(fine)
        coarse_norm = coarse_form.mass_norm_squared(coarse)
        fine_norm = fine_form.mass_norm_squared(fine)
        probes.append({
            "probe": int(index),
            "coarse_energy": coarse_energy,
            "fine_recovery_energy": fine_energy,
            "energy_absolute_delta": abs(fine_energy - coarse_energy),
            "energy_ratio": fine_energy / coarse_energy if coarse_energy > 0 else None,
            "coarse_mass_norm_squared": coarse_norm,
            "fine_mass_norm_squared": fine_norm,
            "mass_norm_absolute_delta": abs(fine_norm - coarse_norm),
        })
    return {
        "status": "FINITE_RECOVERY_PROBES_ONLY",
        "constant_preservation_error": constant_error,
        "probes": probes,
        "scope_boundary": (
            "Finite candidate recovery-sequence diagnostics only; no liminf "
            "proof, dense-core argument, uniform estimate, or Mosco convergence."
        ),
    }
