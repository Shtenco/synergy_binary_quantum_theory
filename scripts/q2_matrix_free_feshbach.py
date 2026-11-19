#!/usr/bin/env python3
"""Matrix-free fibre projection and finite Feshbach operators.

The implementation is graph- and partition-generic. It uses a sparse
combinatorial Laplacian and applies U, U^T, P, and Q by fibre indexing;
it never materializes the fine-by-coarse isometry.
"""
from __future__ import annotations

import numpy as np


class MatrixFreeFeshbach:
    """Matrix-free operators for a graph Laplacian and a fibre partition."""

    def __init__(self, edge_endpoints, fibre_labels, *, validate_edges=True):
        try:
            from scipy.sparse import coo_matrix, diags
        except ImportError as exc:
            raise RuntimeError("matrix-free Feshbach requires scipy") from exc

        edges = np.asarray(edge_endpoints)
        labels = np.asarray(fibre_labels)
        if edges.ndim != 2 or edges.shape[1] != 2:
            raise ValueError("edge_endpoints must have shape (n_edges, 2)")
        if labels.ndim != 1 or labels.size == 0:
            raise ValueError("fibre_labels must be a nonempty one-dimensional array")
        if not np.issubdtype(edges.dtype, np.integer):
            raise ValueError("edge endpoints must be integers")
        if not np.issubdtype(labels.dtype, np.integer):
            raise ValueError("fibre labels must be integers")

        self.n_vertices = int(labels.size)
        self.fibre_labels = labels.astype(np.int64, copy=False)
        if np.any(self.fibre_labels < 0):
            raise ValueError("fibre labels must be non-negative")
        self.n_fibres = int(self.fibre_labels.max()) + 1
        self.fibre_sizes = np.bincount(
            self.fibre_labels, minlength=self.n_fibres
        ).astype(float)
        if np.any(self.fibre_sizes == 0):
            raise ValueError("every coarse fibre must be nonempty")

        edges = edges.astype(np.int64, copy=False)
        if edges.size:
            if np.any(edges < 0) or np.any(edges >= self.n_vertices):
                raise ValueError("edge endpoint outside graph")
            if np.any(edges[:, 0] == edges[:, 1]):
                raise ValueError("self-loops are not supported")
            canonical = np.sort(edges, axis=1)
            if validate_edges and len(np.unique(canonical, axis=0)) != len(edges):
                raise ValueError("duplicate undirected edges are not supported")
            rows = np.concatenate((edges[:, 0], edges[:, 1]))
            cols = np.concatenate((edges[:, 1], edges[:, 0]))
            data = np.ones(rows.size, dtype=float)
            adjacency = coo_matrix(
                (data, (rows, cols)),
                shape=(self.n_vertices, self.n_vertices),
            ).tocsr()
        else:
            adjacency = coo_matrix(
                (self.n_vertices, self.n_vertices), dtype=float
            ).tocsr()

        # QLQ is invertible on range(Q) exactly when ker(L) has no
        # component orthogonal to the fibre-constant subspace. For a graph
        # Laplacian, ker(L) consists of vectors constant on connected
        # components. Thus every connected component must be a union of
        # complete fibres; reject partitions that split a fibre across
        # disconnected components before CG can solve a singular system.
        from scipy.sparse.csgraph import connected_components
        self.graph_component_count, component_labels = connected_components(
            adjacency, directed=False, return_labels=True
        )
        fibre_component_pairs = np.unique(
            np.column_stack((self.fibre_labels, component_labels)), axis=0
        )
        components_per_fibre = np.bincount(
            fibre_component_pairs[:, 0], minlength=self.n_fibres
        )
        if np.any(components_per_fibre > 1):
            raise ValueError(
                "a fibre spans disconnected graph components; QLQ is singular"
            )

        degree = np.asarray(adjacency.sum(axis=1)).ravel()
        self.laplacian = diags(degree, format="csr") - adjacency
        self._inverse_sqrt_sizes = 1.0 / np.sqrt(self.fibre_sizes)
        self._linear_operator = None

    def apply_u(self, coarse):
        coarse = np.asarray(coarse, dtype=float)
        if coarse.shape != (self.n_fibres,):
            raise ValueError("coarse vector has the wrong shape")
        return self._inverse_sqrt_sizes[self.fibre_labels] * coarse[self.fibre_labels]

    def apply_ut(self, fine):
        fine = np.asarray(fine, dtype=float)
        if fine.shape != (self.n_vertices,):
            raise ValueError("fine vector has the wrong shape")
        return np.bincount(
            self.fibre_labels,
            weights=self._inverse_sqrt_sizes[self.fibre_labels] * fine,
            minlength=self.n_fibres,
        )

    def apply_p(self, fine):
        return self.apply_u(self.apply_ut(fine))

    def apply_q(self, fine):
        fine = np.asarray(fine, dtype=float)
        if fine.shape != (self.n_vertices,):
            raise ValueError("fine vector has the wrong shape")
        return fine - self.apply_p(fine)

    def apply_l(self, fine):
        fine = np.asarray(fine, dtype=float)
        if fine.shape != (self.n_vertices,):
            raise ValueError("fine vector has the wrong shape")
        return np.asarray(self.laplacian @ fine).ravel()

    def apply_b(self, coarse):
        return self.apply_q(self.apply_l(self.apply_u(coarse)))

    def apply_bt(self, fine):
        return self.apply_ut(self.apply_l(self.apply_q(fine)))

    def apply_app(self, coarse):
        return self.apply_ut(self.apply_l(self.apply_u(coarse)))

    def apply_q_plus_p(self, fine):
        q_fine = self.apply_q(fine)
        return self.apply_q(self.apply_l(q_fine)) + self.apply_p(fine)

    def _operator(self):
        if self._linear_operator is None:
            from scipy.sparse.linalg import LinearOperator
            self._linear_operator = LinearOperator(
                (self.n_vertices, self.n_vertices),
                matvec=self.apply_q_plus_p,
                rmatvec=self.apply_q_plus_p,
                dtype=float,
            )
        return self._linear_operator

    def solve_q(self, rhs, *, rtol=1e-10, atol=1e-12, maxiter=20000):
        from scipy.sparse.linalg import cg

        rhs = np.asarray(rhs, dtype=float)
        if rhs.shape != (self.n_vertices,):
            raise ValueError("right-hand side has the wrong shape")
        solution, info = cg(
            self._operator(), rhs, rtol=rtol, atol=atol, maxiter=maxiter
        )
        residual = float(np.linalg.norm(self.apply_q_plus_p(solution) - rhs))
        if info != 0:
            raise RuntimeError(f"projected CG did not converge: info={info}")
        if not np.isfinite(residual) or residual > max(atol * 10, rtol * np.linalg.norm(rhs) * 10):
            raise RuntimeError(f"projected CG residual too large: {residual}")
        return solution, {"info": int(info), "residual_norm": residual}

    def apply_feshbach_correction(self, coarse, **solve_options):
        """Apply B^T (QLQ|_Q)^-1 B and return solve diagnostics."""
        fine_rhs = self.apply_b(coarse)
        solution, diagnostics = self.solve_q(fine_rhs, **solve_options)
        return self.apply_bt(solution), diagnostics

    def apply_self_energy_slope(self, coarse, **solve_options):
        """Apply B^T (QLQ|_Q)^-2 B using two projected solves."""
        fine_rhs = self.apply_b(coarse)
        first, first_info = self.solve_q(fine_rhs, **solve_options)
        second, second_info = self.solve_q(first, **solve_options)
        return self.apply_bt(second), {
            "first_solve": first_info,
            "second_solve": second_info,
        }
