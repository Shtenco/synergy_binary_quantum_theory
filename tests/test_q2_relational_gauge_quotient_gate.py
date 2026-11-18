import numpy as np
import pytest

from scripts.q2_relational_graph_qfi_gate import (
    full_half_link_qfi,
    graph_fixtures,
    random_intertwiner,
)
from scripts.q2_relational_gauge_quotient_gate import (
    quotient_diagnostics,
    run,
    two_sided_edge_parameter_map,
    vertex_gauge_generators,
)


def test_global_vertex_gauge_generators_have_three_modes_per_node():
    gauge = vertex_gauge_generators(3)
    assert gauge.shape == (36, 9)
    assert np.allclose(gauge.T @ gauge, 4.0 * np.eye(9))
    assert np.linalg.matrix_rank(gauge) == 9


def test_two_sided_edge_map_keeps_both_endpoint_frames():
    node_count, edges = graph_fixtures()["path_3"]
    mapping = two_sided_edge_parameter_map(node_count, edges)
    assert mapping.shape == (12 * node_count, 6 * len(edges))
    # Each endpoint/axis parameter changes exactly one half-link coordinate.
    assert np.all(np.sum(mapping, axis=0) == 1.0)
    assert np.linalg.matrix_rank(mapping) == 6 * len(edges)


def test_quotient_annihilates_gauss_modes_and_is_section_independent():
    rng = np.random.default_rng(20261118)
    for name in ("path_3", "cycle_4"):
        node_count, edges = graph_fixtures()[name]
        states = [random_intertwiner(rng) for _ in range(node_count)]
        _p, qfi, _sym, _boundary = full_half_link_qfi(states, edges)
        result = quotient_diagnostics(qfi, node_count, reference_legs=(0, 3))
        assert result["gauge_annihilation_residual_frobenius"] < 1e-9
        assert result["quotient_dimension"] == 9 * node_count
        assert result["quotient_is_psd_up_to_roundoff"] is True
        assert result["reference_section_spectra_agree"] is True
        assert result["quotient_nullity_beyond_gauge"] >= 0


def test_reference_leg_must_be_valid():
    with pytest.raises(ValueError, match="reference_leg"):
        from scripts.q2_relational_gauge_quotient_gate import _section_map
        _section_map(2, 4)


def test_end_to_end_path_and_cycle_gate_passes():
    result = run(seed=20261118, trials=3)
    assert result["passed"] is True
    for row in result["graphs"].values():
        assert row["checks"]["global_su2v_gauge_annihilated"]
        assert row["checks"]["reference_section_independence"]
        assert row["checks"]["two_sided_endpoint_parameters_present"]
        assert row["quotient_numerical_rank_min"] <= row["quotient_dimension"]

def test_extra_quotient_nulls_match_projective_boundary_state_stabilizers():
    from scripts.q2_relational_gauge_quotient_gate import (
        projective_tangent_jacobian,
        stabilizer_diagnostics,
    )
    rng = np.random.default_rng(20261118)
    for name in ("path_3", "cycle_4"):
        node_count, edges = graph_fixtures()[name]
        states = [random_intertwiner(rng) for _ in range(node_count)]
        _p, qfi, _sym, _boundary = full_half_link_qfi(states, edges)
        tangent, _p2 = projective_tangent_jacobian(states, edges)
        result = stabilizer_diagnostics(tangent, node_count, qfi)
        assert result["tangent_qfi_rank_agrees"] is True
        assert result["quotient_tangent_nullity"] == 9 * node_count - result["quotient_tangent_rank"]
