#!/usr/bin/env python3
"""Candidate global correspondence: PL tetrahedra -> active Q2 causal edges.

The map is many-to-one. A PL tetrahedron is assigned to one active causal
edge; every PL-dual adjacency must map to two causal edges sharing an endpoint.
Route labels are selected by the colour of the parent face opposite the first
vertex in a barycentric flag. The L/R child is selected by the side assigned
to the parent face opposite the final vertex.

For a shared parent face:
* if its two parent tetrahedra map to distinct causal edges, those edges must
  share exactly one endpoint; each side is forced by which endpoint is shared;
* if they map to the same causal edge, the oriented PL boundary fixes
  opposite sides on the two incident tetrahedra.

This is a finite combinatorial correspondence gate, not a metric, QFI, or
continuum-limit derivation. It is audited only through PL generation 2.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from q2_pl_graph_export import (
    barycentric_subdivision,
    cross_polytope_boundary_4,
    dual_graph_payload,
)
from q2_causal_graph_export import ROUTE_LABELS


def _causal_edges_share_exactly_one_endpoint(edge_a, edge_b):
    """Return whether two distinct endpoint pairs meet at exactly one vertex."""
    return len(set(edge_a) & set(edge_b)) == 1


def _causal_edge_labels_are_incident(edge_id_a, edge_id_b, causal_edges):
    """Same coarse edge is valid; distinct coarse edges must meet at one endpoint."""
    if edge_id_a == edge_id_b:
        return True
    return _causal_edges_share_exactly_one_endpoint(
        causal_edges[edge_id_a], causal_edges[edge_id_b]
    )


def _permutation_data(chain, faces, parent):
    """Recover the ordered parent-vertex permutation from a barycentric flag."""
    permutation = []
    previous = set()
    for face_id in chain:
        current = set(faces[face_id])
        added = current - previous
        if len(added) != 1:
            raise ValueError("barycentric flag is not a maximal simplex chain")
        permutation.append(next(iter(added)))
        previous = current
    if len(permutation) != 4 or set(permutation) != set(parent):
        raise ValueError("flag does not enumerate its parent tetrahedron")
    return permutation


def _vertex_face_colours(dual):
    colours = [dict() for _ in dual["vertices"]]
    for edge in dual["edges"]:
        face = set(edge["shared_face"])
        for endpoint in (edge["source"], edge["target"]):
            tetrahedron = set(dual["vertices"][endpoint]["tetrahedron"])
            omitted = tetrahedron - face
            if len(omitted) != 1:
                raise ValueError("dual edge does not identify one opposite vertex")
            vertex = next(iter(omitted))
            old = colours[endpoint].get(vertex)
            if old is not None and old != edge["route_colour"]:
                raise ValueError("inconsistent route colour on a tetrahedral face")
            colours[endpoint][vertex] = edge["route_colour"]
    if any(sorted(mapping.values()) != [0, 1, 2, 3] for mapping in colours):
        raise ValueError("each tetrahedron must have one face of every route colour")
    return colours


def _boundary_side(tetrahedron, face, orientation):
    omitted = next(iter(set(tetrahedron) - set(face)))
    position = tuple(sorted(tetrahedron)).index(omitted)
    sign = orientation * (-1 if position % 2 else 1)
    return "L" if sign > 0 else "R"


def _parent_face_sides(tetrahedra, parent_map, causal_edges, orientations, dual):
    """Derive L/R for every oriented parent face; return diagnostics too."""
    sides = {}
    same_edge_faces = 0
    distinct_edge_faces = 0
    for edge in dual["edges"]:
        a, b = edge["source"], edge["target"]
        ta = tuple(tetrahedra[a])
        tb = tuple(tetrahedra[b])
        face = tuple(sorted(edge["shared_face"]))
        ea, eb = parent_map[ta], parent_map[tb]
        if ea == eb:
            same_edge_faces += 1
            side_a = _boundary_side(ta, face, orientations[ta])
            side_b = _boundary_side(tb, face, orientations[tb])
            if side_a == side_b:
                raise ValueError("oriented PL boundary sides do not cancel")
        else:
            distinct_edge_faces += 1
            ua, va = causal_edges[ea]
            ub, vb = causal_edges[eb]
            common = set((ua, va)) & set((ub, vb))
            if len(common) != 1:
                raise ValueError("adjacent PL parents map to non-incident causal edges")
            shared = next(iter(common))
            side_a = "L" if ua == shared else "R"
            side_b = "L" if ub == shared else "R"
        sides[(ta, face)] = side_a
        sides[(tb, face)] = side_b
    valid_values = all(side in {"L", "R"} for side in sides.values())
    expected_assignments = 2 * len(dual["edges"])
    return sides, {
        "same_parent_edge_faces": same_edge_faces,
        "distinct_incident_parent_edge_faces": distinct_edge_faces,
        "face_side_assignments": len(sides),
        "expected_face_side_assignments": expected_assignments,
        "side_assignment_values_valid": valid_values,
        "side_assignments_valid": valid_values and len(sides) == expected_assignments,
    }


def _causal_children(causal_edges, generation, next_vertex):
    children = {}
    parent_by_child = {}
    for parent_id, (source, target) in causal_edges.items():
        for label in ROUTE_LABELS:
            midpoint = next_vertex
            next_vertex += 1
            for side, endpoints in (
                ("L", (source, midpoint)),
                ("R", (midpoint, target)),
            ):
                child_id = f"g{generation}:{parent_id}:{label}:{side}"
                children[child_id] = endpoints
                parent_by_child[child_id] = parent_id
    return children, parent_by_child, next_vertex


def _lumpability_diagnostics(tetrahedra, parent_map, dual):
    """Check whether the PL partition is an exact Markov lumping."""
    neighbours = [[] for _ in tetrahedra]
    for edge in dual["edges"]:
        a, b = edge["source"], edge["target"]
        neighbours[a].append(b)
        neighbours[b].append(a)
    fibres = {}
    for index, tet in enumerate(tetrahedra):
        fibres.setdefault(parent_map[tuple(tet)], []).append(index)

    violations = 0
    max_variation = 0
    for source_edge, members in fibres.items():
        target_edges = set()
        counts_by_vertex = []
        for vertex in members:
            counts = Counter(
                parent_map[tuple(tetrahedra[nbr])]
                for nbr in neighbours[vertex]
            )
            counts_by_vertex.append(counts)
            target_edges.update(counts)
        for target_edge in target_edges:
            values = [counts.get(target_edge, 0) for counts in counts_by_vertex]
            variation = max(values) - min(values)
            if variation:
                violations += 1
                max_variation = max(max_variation, variation)
    return {
        "equitable_partition": violations == 0,
        "non_equitable_source_target_pairs": violations,
        "maximum_neighbour_count_variation": max_variation,
    }


def _quotient_operator_diagnostics(generation, tetrahedra, parent_map, causal_edges, dual):
    """Build the sparse Galerkin conductances induced by PL dual adjacency."""
    weights = Counter()
    internal = 0
    for edge in dual["edges"]:
        a = tuple(dual["vertices"][edge["source"]]["tetrahedron"])
        b = tuple(dual["vertices"][edge["target"]]["tetrahedron"])
        ea, eb = parent_map[a], parent_map[b]
        if ea == eb:
            internal += 1
        else:
            weights[tuple(sorted((ea, eb)))] += 1

    nodes = sorted(causal_edges)
    adjacency = {node: set() for node in nodes}
    weighted_degree = {node: 0 for node in nodes}
    sparse_edges = []
    for (a, b), weight in sorted(weights.items()):
        adjacency[a].add(b)
        adjacency[b].add(a)
        weighted_degree[a] += weight
        weighted_degree[b] += weight
        sparse_edges.append({
            "source": a,
            "target": b,
            "integer_multiplicity": weight,
            "normalized_conductance": weight / (24 ** generation),
        })

    seen = set()
    components = 0
    for root in nodes:
        if root in seen:
            continue
        components += 1
        seen.add(root)
        stack = [root]
        while stack:
            current = stack.pop()
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)

    return {
        "quotient_connected": components == 1,
        "quotient_components": components,
        "quotient_weighted_edge_pairs": len(sparse_edges),
        "quotient_cross_fibre_adjacencies": sum(weights.values()),
        "internal_pl_dual_adjacencies": internal,
        "quotient_weighted_degree_min": min(weighted_degree.values()),
        "quotient_weighted_degree_max": max(weighted_degree.values()),
        "quotient_conductances": sparse_edges,
        "conductance_normalization": f"integer multiplicity / 24^{generation}",
        "conductance_semantics": (
            "number of PL-dual edges between the two causal-edge fibres; "
            "Galerkin stiffness weights, not QFI conductances"
        ),
    }


def _jacobi_eigenvalues_symmetric(matrix, tolerance=1e-12, max_sweeps=100):
    """Eigenvalues of a small real symmetric matrix; stdlib-only Jacobi method."""
    n = len(matrix)
    a = [list(map(float, row)) for row in matrix]
    if any(len(row) != n for row in a):
        raise ValueError("symmetric eigensolver requires a square matrix")
    for i in range(n):
        for j in range(i + 1, n):
            if abs(a[i][j] - a[j][i]) > 1e-10:
                raise ValueError("matrix is not symmetric")
            value = 0.5 * (a[i][j] + a[j][i])
            a[i][j] = a[j][i] = value
    for _ in range(max_sweeps):
        changed = False
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) <= tolerance:
                    continue
                changed = True
                tau = (a[q][q] - a[p][p]) / (2.0 * apq)
                sign = 1.0 if tau >= 0 else -1.0
                t = sign / (abs(tau) + (1.0 + tau * tau) ** 0.5)
                cosine = 1.0 / (1.0 + t * t) ** 0.5
                sine = t * cosine
                app, aqq = a[p][p], a[q][q]
                a[p][p] = app - t * apq
                a[q][q] = aqq + t * apq
                a[p][q] = a[q][p] = 0.0
                for k in range(n):
                    if k == p or k == q:
                        continue
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = a[p][k] = cosine * akp - sine * akq
                    a[k][q] = a[q][k] = sine * akp + cosine * akq
        if not changed:
            break
    else:
        residual = max(
            (abs(a[i][j]) for i in range(n) for j in range(i + 1, n)),
            default=0.0,
        )
        if residual > 1e-9:
            raise RuntimeError(f"Jacobi eigensolver did not converge: {residual}")
    return sorted(a[i][i] for i in range(n))


def _symmetric_eigenvalues(matrix):
    """Use LAPACK for large dense matrices, retaining the stdlib Jacobi fallback."""
    n = len(matrix)
    if n >= 128:
        try:
            import numpy as np
        except ImportError:
            return _jacobi_eigenvalues_symmetric(matrix)
        array = np.asarray(matrix, dtype=float)
        if array.ndim != 2 or array.shape != (n, n):
            raise ValueError("symmetric eigensolver requires a square matrix")
        if not np.allclose(array, array.T, rtol=0.0, atol=1e-10):
            raise ValueError("matrix is not symmetric")
        return np.linalg.eigvalsh(array).tolist()
    return _jacobi_eigenvalues_symmetric(matrix)


def _spectral_diagnostics(generation, causal_edges, mass_by_edge, conductances):
    """Spectrum and finite-time heat trace of M^{-1}K, symmetrized by M."""
    nodes = sorted(causal_edges)
    n = len(nodes)
    index = {node: i for i, node in enumerate(nodes)}
    matrix = [[0.0] * n for _ in range(n)]
    scale = 24 ** generation
    for edge in conductances:
        i, j = index[edge["source"]], index[edge["target"]]
        weight = edge["integer_multiplicity"] / scale
        mi, mj = mass_by_edge[edge["source"]], mass_by_edge[edge["target"]]
        if mi <= 0 or mj <= 0:
            raise ValueError("pushforward measure must be strictly positive")
        matrix[i][i] += weight / mi
        matrix[j][j] += weight / mj
        off = -weight / (mi * mj) ** 0.5
        matrix[i][j] += off
        matrix[j][i] += off

    eigenvalues = _symmetric_eigenvalues(matrix)
    zero_count = sum(abs(value) <= 1e-9 for value in eigenvalues)
    positive = [value for value in eigenvalues if value > 1e-9]
    times = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
    heat = []
    for time in times:
        trace = sum(math.exp(-value * time) for value in eigenvalues) / n
        heat.append({"time": time, "normalized_heat_trace": trace})
    dimension = []
    for i in range(1, len(heat) - 1):
        left, right = heat[i - 1], heat[i + 1]
        slope = (
            math.log(right["normalized_heat_trace"])
            - math.log(left["normalized_heat_trace"])
        ) / (math.log(right["time"]) - math.log(left["time"]))
        dimension.append({
            "time": heat[i]["time"],
            "effective_spectral_dimension": -2.0 * slope,
        })
    return {
        "operator_spectrum": eigenvalues,
        "operator_zero_mode_count": zero_count,
        "operator_spectral_gap": min(positive) if positive else None,
        "operator_largest_eigenvalue": max(eigenvalues) if eigenvalues else 0.0,
        "operator_trace": sum(eigenvalues),
        "heat_trace": heat,
        "effective_spectral_dimension": dimension,
        "spectral_interpretation": (
            "finite generalized Galerkin operator spectrum only; a plateau "
            "across refinement generations and a controlled continuum limit "
            "are not established"
        ),
    }


def _certified_q_gap_lower_bound(n, dual):
    """Certified lower bound on lambda_min(QLQ|Q) from graph connectivity.

    For a connected n-vertex graph of diameter D, pathwise Cauchy-Schwarz
    gives lambda_2(L) >= 2 / ((n-1) D). A BFS eccentricity e from any root
    obeys D <= 2e, hence lambda_2(L) >= 1 / ((n-1)e). Since constants lie
    in the fibre-constant P sector, Q is orthogonal to constants and
    lambda_min(QLQ|Q) >= lambda_2(L).
    """
    if n <= 1:
        return None, None
    adjacency = [[] for _ in range(n)]
    for edge in dual["edges"]:
        a, b = edge["source"], edge["target"]
        adjacency[a].append(b)
        adjacency[b].append(a)
    distances = [-1] * n
    distances[0] = 0
    queue = [0]
    for vertex in queue:
        next_distance = distances[vertex] + 1
        for neighbour in adjacency[vertex]:
            if distances[neighbour] < 0:
                distances[neighbour] = next_distance
                queue.append(neighbour)
    if len(queue) != n:
        return None, None
    eccentricity = max(distances)
    if eccentricity <= 0:
        return None, None
    return 1.0 / ((n - 1) * eccentricity), eccentricity


def _feshbach_diagnostics(generation, tetrahedra, parent_map, causal_edges, dual):
    """Compute finite Feshbach data without dense fine-by-coarse embeddings."""
    try:
        import numpy as np
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import LinearOperator, cg, eigsh
        from scipy.linalg import eigh
    except ImportError as exc:
        raise RuntimeError("Feshbach diagnostics require numpy and scipy") from exc

    n = len(tetrahedra)
    certified_gap_lower_bound, root_eccentricity_upper_bound = (
        _certified_q_gap_lower_bound(n, dual)
    )
    nodes = sorted(causal_edges)
    m = len(nodes)
    node_index = {node: index for index, node in enumerate(nodes)}
    fibre_labels = np.fromiter(
        (node_index[parent_map[tuple(tet)]] for tet in tetrahedra),
        dtype=np.int64, count=n,
    )
    fibre_sizes = np.bincount(fibre_labels, minlength=m).astype(float)
    if np.any(fibre_sizes <= 0):
        raise ValueError("every causal coarse fibre must be nonempty")
    inverse_sqrt_sizes = 1.0 / np.sqrt(fibre_sizes)

    # U and U^T are applied by fibre indexing/reduction. Never allocate the
    # n-by-m dense isometry: at generation 3 it would exceed 900 MB.
    def apply_u(coarse):
        return inverse_sqrt_sizes[fibre_labels] * coarse[fibre_labels]

    def apply_ut(fine):
        return np.bincount(
            fibre_labels,
            weights=inverse_sqrt_sizes[fibre_labels] * fine,
            minlength=m,
        )

    def project_q(fine):
        return fine - apply_u(apply_ut(fine))

    rows, cols = [], []
    for edge in dual["edges"]:
        a, b = edge["source"], edge["target"]
        rows.extend((a, b))
        cols.extend((b, a))
    adjacency = csr_matrix(
        (np.ones(len(rows), dtype=float), (rows, cols)), shape=(n, n)
    )

    def laplacian_apply(vector):
        return 4.0 * vector - adjacency @ vector

    A_PP = np.zeros((m, m), dtype=float)
    gram = np.zeros((m, m), dtype=float)
    basis = np.zeros(m, dtype=float)
    for column in range(m):
        basis[column] = 1.0
        u_column = apply_u(basis)
        basis[column] = 0.0
        lu_column = laplacian_apply(u_column)
        b_column = project_q(lu_column)
        A_PP[:, column] = apply_ut(lu_column)
        # B^T B = U^T L B because B=Q L U and B is Q-valued.
        gram[:, column] = apply_ut(laplacian_apply(b_column))

    gram = 0.5 * (gram + gram.T)
    A_PP = 0.5 * (A_PP + A_PP.T)
    leakage_eigenvalues = np.linalg.eigvalsh(gram)
    leakage_frobenius = math.sqrt(max(0.0, float(np.trace(gram))))
    leakage_operator = math.sqrt(max(0.0, float(leakage_eigenvalues[-1])))

    q_sector_gap_estimate = None
    gap_ritz_residual = None
    largest_ritz_residual = None
    q_plus_p_largest_eigenvalue = None
    solve_residual_max = 0.0
    second_solve_residual_max = 0.0
    solve_count = 0
    correction = np.zeros((m, m), dtype=float)
    self_energy_slope = np.zeros((m, m), dtype=float)

    if leakage_frobenius > 1e-14:
        def qlq_plus_p(vector):
            p_vector = apply_u(apply_ut(vector))
            q_vector = vector - p_vector
            return project_q(laplacian_apply(q_vector)) + p_vector

        operator = LinearOperator((n, n), matvec=qlq_plus_p, dtype=float)
        smallest_values, smallest_vectors = eigsh(
            operator, k=1, which="SA", return_eigenvectors=True, tol=1e-12,
            maxiter=200000,
        )
        smallest = float(smallest_values[0])
        gap_ritz_residual = float(np.linalg.norm(
            qlq_plus_p(smallest_vectors[:, 0])
            - smallest * smallest_vectors[:, 0]
        ))
        largest_values, largest_vectors = eigsh(
            operator, k=1, which="LA", return_eigenvectors=True, tol=1e-12,
            maxiter=200000,
        )
        largest = float(largest_values[0])
        largest_ritz_residual = float(np.linalg.norm(
            qlq_plus_p(largest_vectors[:, 0])
            - largest * largest_vectors[:, 0]
        ))
        # P contributes an artificial eigenvalue 1. At the audited levels
        # the lowest eigenvalue is below 1 and therefore belongs to Q.
        if smallest < 1.0 - 1e-8:
            q_sector_gap_estimate = smallest
        q_plus_p_largest_eigenvalue = largest

        for column in range(m):
            basis[column] = 1.0
            u_column = apply_u(basis)
            basis[column] = 0.0
            b_column = project_q(laplacian_apply(u_column))

            solution, info = cg(
                operator, b_column, rtol=1e-10, atol=1e-12, maxiter=20000
            )
            if info != 0:
                raise RuntimeError(
                    f"Feshbach CG failed at generation {generation}, "
                    f"column {column}, info={info}"
                )
            solve_residual_max = max(
                solve_residual_max,
                float(np.linalg.norm(qlq_plus_p(solution) - b_column)),
            )
            correction[:, column] = apply_ut(
                laplacian_apply(project_q(solution))
            )

            second_solution, info = cg(
                operator, solution, rtol=1e-10, atol=1e-12, maxiter=20000
            )
            if info != 0:
                raise RuntimeError(
                    f"Feshbach second-resolvent CG failed at generation "
                    f"{generation}, column {column}, info={info}"
                )
            second_solve_residual_max = max(
                second_solve_residual_max,
                float(np.linalg.norm(qlq_plus_p(second_solution) - solution)),
            )
            self_energy_slope[:, column] = apply_ut(
                laplacian_apply(project_q(second_solution))
            )
            solve_count += 1

        correction = 0.5 * (correction + correction.T)
        self_energy_slope = 0.5 * (
            self_energy_slope + self_energy_slope.T
        )

    correction_eigenvalues = np.linalg.eigvalsh(correction)
    slope_eigenvalues = np.linalg.eigvalsh(self_energy_slope)
    effective = A_PP - correction
    effective = 0.5 * (effective + effective.T)
    effective_eigenvalues = np.linalg.eigvalsh(effective)
    # Linearize the energy-dependent Schur complement:
    # A_eff(0) v = lambda [I + Sigma'(0)] v + O(lambda^2).
    effective_metric = np.eye(m) + self_energy_slope
    effective_metric = 0.5 * (effective_metric + effective_metric.T)
    metric_eigenvalues = np.linalg.eigvalsh(effective_metric)
    generalized_eigenvalues = eigh(
        effective, effective_metric, eigvals_only=True, check_finite=True
    )
    app_norm = float(np.linalg.norm(A_PP, ord=2))
    correction_norm = float(np.linalg.norm(correction, ord=2))
    return {
        "method": "matrix-free fibre U/U^T; sparse PL adjacency; CG on QLQ+P",
        "pl_vertices": n,
        "coarse_modes": m,
        "fine_by_coarse_dense_embedding_allocated": False,
        "leakage_frobenius_norm": leakage_frobenius,
        "leakage_operator_norm": leakage_operator,
        "leakage_rms_per_coarse_mode": leakage_frobenius / math.sqrt(m),
        "A_PP_operator_norm": app_norm,
        "leakage_to_galerkin_operator_ratio": (
            leakage_operator / app_norm if app_norm > 0 else None
        ),
        "feshbach_correction_operator_norm": correction_norm,
        "feshbach_correction_frobenius_norm": float(np.linalg.norm(correction, ord="fro")),
        "feshbach_to_galerkin_operator_ratio": (
            correction_norm / app_norm if app_norm > 0 else None
        ),
        "self_energy_slope_operator_norm": float(np.linalg.norm(self_energy_slope, ord=2)),
        "self_energy_slope_frobenius_norm": float(np.linalg.norm(self_energy_slope, ord="fro")),
        "self_energy_slope_trace": float(np.trace(self_energy_slope)),
        "self_energy_slope_min_eigenvalue": float(slope_eigenvalues[0]),
        "self_energy_slope_max_eigenvalue": float(slope_eigenvalues[-1]),
        "self_energy_slope_to_identity_ratio": float(np.linalg.norm(self_energy_slope, ord=2)),
        "max_second_resolvent_solve_residual": second_solve_residual_max,
        "feshbach_correction_min_eigenvalue": float(correction_eigenvalues[0]),
        "feshbach_correction_max_eigenvalue": float(correction_eigenvalues[-1]),
        "effective_operator_min_eigenvalue": float(effective_eigenvalues[0]),
        "effective_operator_gap": (
            float(effective_eigenvalues[1]) if len(effective_eigenvalues) > 1 else None
        ),
        "effective_operator_zero_mode_count": int(
            np.sum(np.abs(effective_eigenvalues) <= 1e-9)
        ),
        "effective_operator_max_eigenvalue": float(effective_eigenvalues[-1]),
        "effective_metric_min_eigenvalue": float(metric_eigenvalues[0]),
        "effective_metric_max_eigenvalue": float(metric_eigenvalues[-1]),
        "linearized_generalized_spectrum": generalized_eigenvalues.tolist(),
        "linearized_generalized_gap": (
            float(generalized_eigenvalues[1])
            if len(generalized_eigenvalues) > 1 else None
        ),
        "linearized_generalized_largest_eigenvalue": float(generalized_eigenvalues[-1]),
        "q_sector_gap_estimate": q_sector_gap_estimate,
        "q_sector_gap_ritz_residual": gap_ritz_residual,
        "q_sector_gap_certified_lower_bound": certified_gap_lower_bound,
        "root_eccentricity_upper_bound": root_eccentricity_upper_bound,
        "q_plus_p_largest_eigenvalue": q_plus_p_largest_eigenvalue,
        "q_plus_p_largest_ritz_residual": largest_ritz_residual,
        "q_plus_p_condition_estimate": (
            q_plus_p_largest_eigenvalue / q_sector_gap_estimate
            if q_sector_gap_estimate is not None else None
        ),
        "self_energy_linearization_remainder_coefficient_using_gap_estimate": (
            leakage_operator ** 2 / (q_sector_gap_estimate ** 3)
            if q_sector_gap_estimate is not None else None
        ),
        "self_energy_linearization_remainder_coefficient_certified_upper_bound": (
            leakage_operator ** 2 / (certified_gap_lower_bound ** 3)
            if certified_gap_lower_bound is not None else None
        ),
        "remainder_bound_certified": certified_gap_lower_bound is not None,
        "max_linear_solve_residual": solve_residual_max,
        "linear_solve_columns": solve_count,
        "all_linear_solves_converged": True,
        "interpretation": (
            "finite-level zero-energy Feshbach correction for the unweighted "
            "PL dual Laplacian; not a continuum or physical QFI result"
        ),
    }


def _summarise(generation, tetrahedra, parent_map, causal_edges, dual, compute_feshbach=False):
    fibres = Counter(parent_map[tuple(tet)] for tet in tetrahedra)
    adjacency_violations = []
    for edge in dual["edges"]:
        a = tuple(dual["vertices"][edge["source"]]["tetrahedron"])
        b = tuple(dual["vertices"][edge["target"]]["tetrahedron"])
        ea, eb = parent_map[a], parent_map[b]
        endpoints_a, endpoints_b = causal_edges[ea], causal_edges[eb]
        if not _causal_edge_labels_are_incident(ea, eb, causal_edges):
            adjacency_violations.append({
                "pl_dual_source": edge["source"],
                "pl_dual_target": edge["target"],
                "causal_source_edge": ea,
                "causal_target_edge": eb,
                "causal_source_endpoints": list(endpoints_a),
                "causal_target_endpoints": list(endpoints_b),
                "shared_endpoint_count": len(set(endpoints_a) & set(endpoints_b)),
            })
    adjacency_ok = not adjacency_violations
    histogram = Counter(fibres.values())
    measure_scale = 24 ** generation
    lumpability = _lumpability_diagnostics(tetrahedra, parent_map, dual)
    quotient = _quotient_operator_diagnostics(
        generation, tetrahedra, parent_map, causal_edges, dual
    )
    mass_by_edge = {
        edge: count / (24 ** generation) for edge, count in fibres.items()
    }
    spectrum = _spectral_diagnostics(
        generation, causal_edges, mass_by_edge, quotient["quotient_conductances"]
    )
    return {
        "generation": generation,
        "pl_tetrahedra": len(tetrahedra),
        "causal_active_edges": len(causal_edges),
        "surjective": set(fibres) == set(causal_edges),
        "fibre_min": min(fibres.values()),
        "fibre_max": max(fibres.values()),
        "fibre_mean": len(tetrahedra) / len(causal_edges),
        "fibre_size_histogram": {str(k): histogram[k] for k in sorted(histogram)},
        "pushforward_measure_total": sum(value / measure_scale for value in fibres.values()),
        "pushforward_measure_min": min(value / measure_scale for value in fibres.values()),
        "pushforward_measure_max": max(value / measure_scale for value in fibres.values()),
        "pl_dual_adjacencies": len(dual["edges"]),
        "adjacency_preserved": adjacency_ok,
        "adjacency_violation_count": len(adjacency_violations),
        "adjacency_violation_examples": adjacency_violations[:10],
        **lumpability,
        **quotient,
        **spectrum,
        **({"feshbach": _feshbach_diagnostics(
            generation, tetrahedra, parent_map, causal_edges, dual
        )} if compute_feshbach else {}),
    }

def run(max_generation=2, compute_feshbach=False):
    if not 0 <= max_generation <= 2:
        raise ValueError("audited PL exporter supports generations 0..2")

    tetrahedra = cross_polytope_boundary_4()
    parent_map = {tuple(tet): "e0" for tet in tetrahedra}
    orientations = {
        tuple(tet): (-1 if sum(vertex % 2 for vertex in tet) % 2 else 1)
        for tet in tetrahedra
    }
    causal_edges = {"e0": (0, 1)}
    next_vertex = 2
    previous_faces = None
    chain_by_tet = None
    levels = []
    refinement_checks = []
    face_diagnostics = []
    route_colour_diagnostics = []
    measure_checks = []

    for generation in range(max_generation + 1):
        dual = dual_graph_payload(
            generation, tetrahedra, chain_by_tet, previous_faces
        )
        vertex_colours = _vertex_face_colours(dual)
        route_colours_valid = all(
            sorted(mapping.values()) == [0, 1, 2, 3]
            for mapping in vertex_colours
        )
        route_colour_diagnostics.append({
            "generation": generation,
            "vertices_checked": len(vertex_colours),
            "all_four_route_colours_present_per_vertex": route_colours_valid,
        })
        levels.append(_summarise(
            generation, tetrahedra, parent_map, causal_edges, dual,
            compute_feshbach=compute_feshbach,
        ))
        if generation == max_generation:
            break

        face_sides, face_info = _parent_face_sides(
            tetrahedra, parent_map, causal_edges, orientations, dual
        )
        face_diagnostics.append({
            "generation": generation,
            **face_info,
        })

        children, faces, child_chains = barycentric_subdivision(tetrahedra)
        parent_index = {tuple(tet): i for i, tet in enumerate(tetrahedra)}
        child_causal_edges, parent_by_child, next_vertex = _causal_children(
            causal_edges, generation + 1, next_vertex
        )
        child_map = {}
        child_orientations = {}
        naturality_ok = True

        for child in children:
            chain = child_chains[tuple(child)]
            parent = tuple(faces[chain[-1]])
            pidx = parent_index[parent]
            permutation = _permutation_data(chain, faces, parent)
            route_colour = vertex_colours[pidx][permutation[0]]
            omitted_last = permutation[-1]
            parent_face = tuple(sorted(v for v in parent if v != omitted_last))
            side = face_sides[(parent, parent_face)]
            label = ROUTE_LABELS[route_colour]
            causal_id = f"g{generation + 1}:{parent_map[parent]}:{label}:{side}"
            child_map[tuple(child)] = causal_id
            if causal_id not in child_causal_edges:
                naturality_ok = False

            positions = [tuple(sorted(parent)).index(v) for v in permutation]
            inversions = sum(
                positions[i] > positions[j]
                for i in range(4) for j in range(i + 1, 4)
            )
            child_orientations[tuple(child)] = (
                orientations[parent] * (-1 if inversions % 2 else 1)
            )

        refinement_checks.append({
            "generation": generation + 1,
            "every_child_maps_to_a_child_of_its_parent": naturality_ok,
        })

        # Push forward PL counting measure, normalized by the 24 children
        # created per parent tetrahedron. The mass of each parent causal-edge
        # fibre must equal the sum of the masses of its eight child fibres.
        parent_fibres = Counter(parent_map[tuple(tet)] for tet in tetrahedra)
        child_fibres = Counter(child_map[tuple(tet)] for tet in children)
        child_mass_by_parent = Counter()
        for child_edge, count in child_fibres.items():
            child_mass_by_parent[parent_by_child[child_edge]] += count
        mass_errors = {}
        for parent_edge, parent_count in parent_fibres.items():
            parent_mass = parent_count / (24 ** generation)
            child_mass = child_mass_by_parent[parent_edge] / (24 ** (generation + 1))
            if abs(parent_mass - child_mass) > 1e-12:
                mass_errors[parent_edge] = {
                    "parent_mass": parent_mass,
                    "child_mass": child_mass,
                }
        measure_checks.append({
            "generation": generation + 1,
            "parent_mass_conserved": not mass_errors,
            "parent_fibres_checked": len(parent_fibres),
            "mass_errors": mass_errors,
        })
        # Carry the refined complex forward; its level summary is built once
        # at the start of the next iteration.
        tetrahedra = children
        parent_map = child_map
        orientations = child_orientations
        causal_edges = child_causal_edges
        previous_faces = faces
        chain_by_tet = child_chains

    checks = {
        "all_levels_surjective": all(row["surjective"] for row in levels),
        "all_levels_preserve_dual_adjacency": all(
            row["adjacency_preserved"] for row in levels
        ),
        "all_fibres_nonempty": all(row["fibre_min"] > 0 for row in levels),
        "pushforward_measure_mass_is_16": all(
            abs(row["pushforward_measure_total"] - 16.0) < 1e-12
            for row in levels
        ),
        "parent_pushforward_measure_is_refinement_consistent": all(
            row["parent_mass_conserved"] for row in measure_checks
        ),
        "refinement_naturality": all(
            row["every_child_maps_to_a_child_of_its_parent"]
            for row in refinement_checks
        ),
        "all_parent_face_side_assignments_valid": all(
            row["side_assignments_valid"] for row in face_diagnostics
        ),
        "four_route_face_colours_present": all(
            row["all_four_route_colours_present_per_vertex"]
            for row in route_colour_diagnostics
        ),
    }
    return {
        "status": "candidate q2 causal-to-PL edge correspondence",
        "passed": all(checks.values()),
        "map_semantics": (
            "PL tetrahedra map many-to-one onto active causal edges; every PL "
            "dual adjacency maps to endpoint incidence of causal edges"
        ),
        "side_rule": (
            "distinct parent causal edges: choose the child side containing "
            "their shared endpoint; equal parent edge: use opposite oriented "
            "PL boundary sides"
        ),
        "levels": levels,
        "parent_face_diagnostics": face_diagnostics,
        "route_colour_diagnostics": route_colour_diagnostics,
        "refinement_checks": refinement_checks,
        "measure_checks": measure_checks,
        "measure_definition": (
            "mu_g(e) = |Phi_g^{-1}(e)| / 24^g; total mass is 16 and "
            "parent mass equals the sum over its eight child-edge masses"
        ),
        "checks": checks,
        "claim_boundary": (
            "Finite combinatorial graph correspondence only. Nonuniform fibres "
            "at generation 2 are reported, not hidden or reweighted. This does "
            "not derive a physical continuum measure, QFI conductance, spectral dimension, "
            "semiclassical dynamics, or spacetime."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-generation", type=int, default=2)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--feshbach", action="store_true",
                        help="compute matrix-free finite-level Feshbach corrections")
    args = parser.parse_args()
    result = run(args.max_generation, compute_feshbach=args.feshbach)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
