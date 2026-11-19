import numpy as np
import pytest

from scripts.q2_pl_geometric_mesh_audit import (
    build_geometry,
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
    assert np.allclose(volumes, 1.0 / 3.0)


def test_barycentric_refinement_conserves_pl_volume_and_shrinks_mesh():
    _, tetrahedra, _, levels = build_geometry(2)
    assert [row["tetrahedra"] for row in levels] == [16, 384, 9216]
    assert all(row["pl_volume_total"] == pytest.approx(16.0 / 3.0, abs=1e-10)
               for row in levels)
    assert levels[0]["mesh_diameter_max"] > levels[1]["mesh_diameter_max"]
    assert levels[1]["mesh_diameter_max"] > levels[2]["mesh_diameter_max"]
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
