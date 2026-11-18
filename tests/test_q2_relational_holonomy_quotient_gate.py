import numpy as np

from scripts.q2_relational_holonomy_quotient_gate import (
    adjoint_matrix,
    diagnostics,
    endpoint_to_holonomy_jacobian,
    fixture_holonomies,
    graph_fixtures,
    holonomy_gauge_jacobian,
    loop_holonomy,
    run,
    su2_from_vector,
    transform_holonomies,
    vertex_gauge_embedding,
)


def test_adjoint_representation_is_orthogonal_and_composes():
    g = su2_from_vector(np.array([0.23, -0.31, 0.17]))
    h = su2_from_vector(np.array([-0.11, 0.07, 0.29]))
    assert np.allclose(adjoint_matrix(g).T @ adjoint_matrix(g), np.eye(3), atol=2e-10)
    assert np.allclose(adjoint_matrix(g @ h), adjoint_matrix(g) @ adjoint_matrix(h), atol=2e-10)


def test_six_endpoint_coordinates_map_onto_three_holonomy_tangents():
    for name, (nodes, edges) in graph_fixtures().items():
        holonomies = fixture_holonomies(name)
        mapping = endpoint_to_holonomy_jacobian(holonomies)
        assert mapping.shape == (3 * len(edges), 6 * len(edges))
        assert np.linalg.matrix_rank(mapping, tol=1e-9) == 3 * len(edges)


def test_vertex_gauge_embedding_composes_to_correct_adjoint_incidence_map():
    for name, (nodes, edges) in graph_fixtures().items():
        holonomies = fixture_holonomies(name)
        endpoint = endpoint_to_holonomy_jacobian(holonomies)
        embedding = vertex_gauge_embedding(nodes, edges)
        composed = endpoint @ embedding
        expected = holonomy_gauge_jacobian(nodes, edges, holonomies)
        assert np.allclose(composed, expected, atol=1e-12)
        assert composed.shape == (3 * len(edges), 3 * nodes)


def test_tree_and_single_cycle_have_different_holonomy_moduli():
    path_nodes, path_edges = graph_fixtures()["path_3"]
    cycle_nodes, cycle_edges = graph_fixtures()["cycle_4"]
    path = diagnostics("path_3", path_nodes, path_edges, fixture_holonomies("path_3"))
    cycle = diagnostics("cycle_4", cycle_nodes, cycle_edges, fixture_holonomies("cycle_4"))
    assert path["cycle_rank"] == 0
    assert path["holonomy_mod_vertex_gauge_tangent_dimension"] == 0
    assert path["vertex_gauge_stabilizer_dimension"] == 3
    assert cycle["cycle_rank"] == 1
    assert cycle["holonomy_mod_vertex_gauge_tangent_dimension"] == 1
    assert cycle["vertex_gauge_stabilizer_dimension"] == 1


def test_finite_vertex_gauge_transformation_conjugates_cycle_holonomy():
    nodes, edges = graph_fixtures()["cycle_4"]
    holonomies = fixture_holonomies("cycle_4")
    transforms = [
        su2_from_vector(np.array([0.1, 0.2, -0.1])),
        su2_from_vector(np.array([-0.2, 0.1, 0.07])),
        su2_from_vector(np.array([0.04, -0.13, 0.19])),
        su2_from_vector(np.array([0.08, 0.03, -0.17])),
    ]
    transformed = transform_holonomies(holonomies, edges, transforms)
    before = loop_holonomy(holonomies)
    after = loop_holonomy(transformed)
    expected = transforms[0] @ before @ transforms[0].conj().T
    assert np.allclose(after, expected, atol=2e-10)
    assert np.isclose(np.trace(after), np.trace(before), atol=2e-10)


def test_end_to_end_holonomy_quotient_gate_passes():
    result = run()
    assert result["passed"] is True
    assert all(result["checks"].values())
