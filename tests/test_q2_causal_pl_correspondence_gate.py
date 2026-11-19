from scripts.q2_causal_pl_correspondence_gate import run


def test_q2_causal_to_pl_correspondence_preserves_adjacency_and_refinement():
    result = run(max_generation=2)
    print("ADJACENCY_DIAGNOSTICS", [(row["generation"], row["adjacency_violation_count"], row["adjacency_violation_examples"]) for row in result["levels"]])
    assert result["passed"] is True, {"checks": result["checks"]}
    assert result["checks"]["all_levels_preserve_dual_adjacency"] is True
    assert result["checks"]["refinement_naturality"] is True
    assert [
        (row["pl_tetrahedra"], row["causal_active_edges"])
        for row in result["levels"]
    ] == [(16, 1), (384, 8), (9216, 64)]


def test_q2_causal_to_pl_correspondence_fibre_multiplicities_are_explicit():
    result = run(max_generation=2)
    level0, level1, level2 = result["levels"]
    assert level0["fibre_size_histogram"] == {"16": 1}
    assert level1["fibre_size_histogram"] == {"48": 8}
    assert level2["fibre_size_histogram"] == {
        "64": 8,
        "112": 8,
        "128": 8,
        "144": 16,
        "160": 8,
        "176": 8,
        "224": 8,
    }
    assert level2["fibre_min"] == 64
    assert level2["fibre_max"] == 224
    assert level2["fibre_mean"] == 144


def test_q2_causal_to_pl_correspondence_does_not_claim_physical_weights():
    result = run(max_generation=2)
    assert "does not derive a physical continuum measure" in result["claim_boundary"]



def test_q2_pushforward_measure_is_conserved_under_refinement():
    result = run(max_generation=2)
    assert result["checks"]["pushforward_measure_mass_is_16"] is True
    assert result["checks"]["parent_pushforward_measure_is_refinement_consistent"] is True
    assert [row["pushforward_measure_total"] for row in result["levels"]] == [16.0, 16.0, 16.0]
    assert all(row["parent_mass_conserved"] for row in result["measure_checks"])


def test_q2_causal_quotient_is_not_an_exact_markov_lumping():
    result = run(max_generation=2)
    level0, level1, level2 = result["levels"]
    assert level0["equitable_partition"] is True
    assert level1["equitable_partition"] is False
    assert level1["non_equitable_source_target_pairs"] == 40
    assert level2["equitable_partition"] is False
    assert level2["non_equitable_source_target_pairs"] == 408
    assert level2["maximum_neighbour_count_variation"] == 2



def test_q2_galerkin_quotient_conductances_are_derived_from_pl_adjacency():
    result = run(max_generation=2)
    level0, level1, level2 = result["levels"]
    assert (level0["quotient_weighted_edge_pairs"], level0["internal_pl_dual_adjacencies"]) == (0, 32)
    assert (
        level1["quotient_weighted_edge_pairs"],
        level1["quotient_cross_fibre_adjacencies"],
        level1["internal_pl_dual_adjacencies"],
        level1["quotient_weighted_degree_min"],
        level1["quotient_weighted_degree_max"],
    ) == (16, 512, 256, 128, 128)
    assert (
        level2["quotient_weighted_edge_pairs"],
        level2["quotient_cross_fibre_adjacencies"],
        level2["internal_pl_dual_adjacencies"],
        level2["quotient_weighted_degree_min"],
        level2["quotient_weighted_degree_max"],
    ) == (180, 12032, 6400, 192, 512)
    assert level1["quotient_connected"] is True
    assert level2["quotient_connected"] is True
    assert all(edge["integer_multiplicity"] > 0 for edge in level2["quotient_conductances"])
    assert "not QFI conductances" in level2["conductance_semantics"]
    first = level2["quotient_conductances"][0]
    assert first["normalized_conductance"] == first["integer_multiplicity"] / (24 ** 2)
    assert level2["conductance_normalization"] == "integer multiplicity / 24^2"



def test_q2_galerkin_operator_has_one_zero_mode_and_positive_spectral_gap():
    result = run(max_generation=2)
    for row in result["levels"][1:]:
        assert row["operator_zero_mode_count"] == 1
        assert row["operator_spectral_gap"] is not None
        assert row["operator_spectral_gap"] > 0
        assert row["operator_largest_eigenvalue"] >= row["operator_spectral_gap"]
        assert min(row["operator_spectrum"]) >= -1e-9
        assert len(row["operator_spectrum"]) == row["causal_active_edges"]


def test_q2_finite_heat_trace_is_monotone_but_not_claimed_as_dimension_limit():
    result = run(max_generation=2)
    for row in result["levels"][1:]:
        traces = [item["normalized_heat_trace"] for item in row["heat_trace"]]
        assert all(a >= b for a, b in zip(traces, traces[1:]))
        assert all(item["effective_spectral_dimension"] >= -1e-8
                   for item in row["effective_spectral_dimension"])
        assert "continuum limit are not established" in row["spectral_interpretation"]


def test_q2_jacobi_eigensolver_matches_exact_two_by_two_spectrum():
    from scripts.q2_causal_pl_correspondence_gate import _jacobi_eigenvalues_symmetric

    values = _jacobi_eigenvalues_symmetric([[2.0, 1.0], [1.0, 2.0]])
    assert len(values) == 2
    assert abs(values[0] - 1.0) < 1e-12
    assert abs(values[1] - 3.0) < 1e-12


def test_q2_jacobi_eigensolver_rejects_non_symmetric_input():
    import pytest
    from scripts.q2_causal_pl_correspondence_gate import _jacobi_eigenvalues_symmetric

    with pytest.raises(ValueError, match="matrix is not symmetric"):
        _jacobi_eigenvalues_symmetric([[1.0, 2.0], [0.0, 1.0]])


def test_q2_generation_zero_spectrum_has_only_constant_zero_mode():
    result = run(max_generation=0)
    row = result["levels"][0]
    assert row["operator_spectrum"] == [0.0]
    assert row["operator_zero_mode_count"] == 1
    assert row["operator_spectral_gap"] is None
    assert all(item["normalized_heat_trace"] == 1.0 for item in row["heat_trace"])

def test_q2_large_symmetric_spectrum_uses_lapack_backend():
    import numpy as np
    from scripts.q2_causal_pl_correspondence_gate import _symmetric_eigenvalues

    matrix = np.diag(np.arange(128, dtype=float))
    values = _symmetric_eigenvalues(matrix.tolist())
    assert np.allclose(values, np.arange(128, dtype=float), atol=1e-12)


def test_q2_large_symmetric_eigensolver_rejects_asymmetric_input():
    import numpy as np
    import pytest
    from scripts.q2_causal_pl_correspondence_gate import _symmetric_eigenvalues

    matrix = np.eye(128)
    matrix[0, 1] = 0.5
    with pytest.raises(ValueError, match="matrix is not symmetric"):
        _symmetric_eigenvalues(matrix.tolist())

def test_q2_causal_edge_incidence_requires_exactly_one_shared_endpoint():
    from scripts.q2_causal_pl_correspondence_gate import (
        _causal_edges_share_exactly_one_endpoint,
    )

    assert _causal_edges_share_exactly_one_endpoint((0, 1), (1, 2)) is True
    assert _causal_edges_share_exactly_one_endpoint((0, 1), (2, 3)) is False
    # Distinct edge labels must not hide a duplicate endpoint pair.
    assert _causal_edges_share_exactly_one_endpoint((0, 1), (0, 1)) is False
    assert _causal_edges_share_exactly_one_endpoint((0, 1), (1, 0)) is False


def test_q2_correspondence_summary_gates_are_computed_from_diagnostics():
    result = run(max_generation=2)

    assert result["checks"]["all_parent_face_side_assignments_valid"] is True
    assert result["checks"]["four_route_face_colours_present"] is True
    assert len(result["route_colour_diagnostics"]) == 3
    assert all(
        row["all_four_route_colours_present_per_vertex"]
        and row["vertices_checked"] > 0
        for row in result["route_colour_diagnostics"]
    )
    assert all(
        row["side_assignments_valid"]
        and row["face_side_assignments"] == row["expected_face_side_assignments"]
        for row in result["parent_face_diagnostics"]
    )

