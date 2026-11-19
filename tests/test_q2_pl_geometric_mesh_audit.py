import numpy as np
import pytest

from scripts.q2_pl_geometric_mesh_audit import (
    assess_shape_regularity,
    build_coarse_spectral_sequence,
    build_geometry,
    coarse_finite_volume_operator,
    geometric_dual_weights,
    initial_reference_embedding,
    mesh_metrics,
)


def test_cross_polytope_reference_boundary_has_volume_16_over_3():
    coordinates, tetrahedra = initial_reference_embedding()
    metrics, volumes = mesh_metrics(coordinates, tetrahedra)
    assert metrics["tetrahedra"] == 16
    assert metrics["pl_volume_total"] == pytest.approx(16.0 / 3.0, abs=1e-12)
    assert metrics["cell_volume_min"] == pytest.approx(1.0 / 3.0, abs=1e-12)
    assert metrics["mesh_diameter_max"] == pytest.approx(np.sqrt(2.0), abs=1e-12)
    assert metrics["shape_quality_min"] == pytest.approx(1.0, abs=1e-12)
    assert metrics["max_diameter_to_inradius"] > 0
    assert np.allclose(volumes, 1.0 / 3.0)


def test_barycentric_refinement_conserves_pl_volume_and_shrinks_mesh():
    _, tetrahedra, _, levels = build_geometry(2)
    assert [row["tetrahedra"] for row in levels] == [16, 384, 9216]
    assert all(row["pl_volume_total"] == pytest.approx(16.0 / 3.0, abs=1e-10)
               for row in levels)
    assert levels[0]["mesh_diameter_max"] > levels[1]["mesh_diameter_max"]
    assert levels[1]["mesh_diameter_max"] > levels[2]["mesh_diameter_max"]
    assert all(row["shape_quality_min"] > 0 for row in levels)
    assert levels[2]["shape_quality_min"] <= levels[0]["shape_quality_min"]
    assert tetrahedra.shape == (9216, 4)


def test_generation_zero_geometric_dual_conductances_are_positive():
    from scripts.q2_pl_graph_export import dual_graph_payload

    coordinates, tetrahedra = initial_reference_embedding()
    metrics, volumes = mesh_metrics(coordinates, tetrahedra)
    dual = dual_graph_payload(0, tetrahedra)
    endpoints = np.asarray(
        [(edge["source"], edge["target"]) for edge in dual["edges"]],
        dtype=np.uint32,
    )
    weights, areas, distances = geometric_dual_weights(
        coordinates, tetrahedra, endpoints, volumes, chunk_rows=7
    )
    assert len(endpoints) == 32
    assert np.all(weights > 0)
    assert np.all(areas > 0)
    assert np.all(distances > 0)
    assert np.all(np.isfinite(weights))

def test_stable_child_order_matches_exporter_id_contract():
    import itertools
    from scripts.q2_pl_graph_export import all_nonempty_faces
    from scripts.q2_pl_generation3_stream import _child_key
    from scripts.q2_pl_geometric_mesh_audit import stable_child_order

    parent = [(0, 1, 2, 3)]
    faces = all_nonempty_faces(parent)
    face_id = {face: index for index, face in enumerate(faces)}
    expected = np.asarray([
        _child_key(parent[0], tuple(parent[0][i] for i in permutation), face_id)
        for permutation in itertools.permutations(range(4))
    ], dtype=np.uint32)
    assert np.array_equal(stable_child_order(parent, faces), expected)

def test_shape_regularity_assessment_does_not_hide_quality_degradation():
    levels = [
        {"shape_quality_min": 1.0, "max_diameter_to_inradius": 4.9},
        {"shape_quality_min": 0.6, "max_diameter_to_inradius": 10.0},
        {"shape_quality_min": 0.2, "max_diameter_to_inradius": 30.0},
    ]
    report = assess_shape_regularity(levels)
    assert report["shape_quality_degradation_detected"]
    assert report["aspect_ratio_growth_detected"]
    assert not report["uniform_shape_regularity_established"]
    assert "cannot be assumed" in report["warning"]

def test_coarse_geometric_operator_has_mass_weighted_constant_nullmode():
    edges = np.asarray([[0, 1], [1, 2]], dtype=np.uint32)
    conductances = np.asarray([2.0, 3.0])
    labels = np.asarray([0, 1, 1], dtype=np.uint16)
    volumes = np.asarray([0.25, 0.25, 0.25])
    operator, mass = coarse_finite_volume_operator(
        edges, conductances, labels, volumes
    )
    constant = np.sqrt(mass)
    constant /= np.linalg.norm(constant)
    assert np.linalg.norm(operator @ constant) < 1e-12
    assert np.linalg.eigvalsh(operator)[0] >= -1e-12
    assert mass.sum() == pytest.approx(0.75)


def test_geometric_galerkin_spectrum_is_built_for_g0_to_g2():
    reports = build_coarse_spectral_sequence()
    assert [row["generation"] for row in reports] == [0, 1, 2]
    assert [row["coarse_dimension"] for row in reports] == [1, 8, 64]
    assert all(row["coarse_mass_total"] == pytest.approx(16.0 / 3.0, abs=1e-10)
               for row in reports)
    assert [row["zero_modes_at_relative_tolerance"] for row in reports] == [1, 1, 1]
    assert reports[1]["smallest_positive_eigenvalue"] > 0
    assert reports[2]["smallest_positive_eigenvalue"] > 0
