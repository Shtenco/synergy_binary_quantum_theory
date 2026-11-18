from scripts.q2_pl_graph_export import run


def test_q2_pl_graph_export_matches_frozen_refinement_counts():
    result = run(max_generation=2)
    assert result["passed"] is True
    assert [
        (row["vertex_count"], row["edge_count"])
        for row in result["generations"]
    ] == [(8, 24), (80, 464), (1696, 10912)]
    assert all(row["component_count"] == 1 for row in result["generations"])


def test_q2_pl_graph_export_preserves_canonical_refinement_labels():
    result = run(max_generation=2)
    level1, level2 = result["generations"][1:]
    assert len(level1["refinement_map"]) == 80
    assert len(level2["refinement_map"]) == 1696
    assert all(
        item["child_vertex"] == index
        and 1 <= len(item["parent_simplex"]) <= 4
        for index, item in enumerate(level1["refinement_map"])
    )
    assert all(
        item["child_vertex"] == index
        and 1 <= len(item["parent_simplex"]) <= 4
        for index, item in enumerate(level2["refinement_map"])
    )


def test_q2_pl_graph_export_does_not_fabricate_quantum_data():
    result = run(max_generation=2)
    for row in result["generations"]:
        assert row["quantum_state_assignment"] is None
        assert row["qfi_edge_conductances"] is None
        assert len(row["graph_sha256"]) == 64



def test_q2_pl_export_includes_four_regular_tetrahedral_dual_graphs():
    result = run(max_generation=2)
    assert [
        (row["dual_graph"]["vertex_count"], row["dual_graph"]["edge_count"])
        for row in result["generations"]
    ] == [(16, 32), (384, 768), (9216, 18432)]
    for row in result["generations"]:
        dual = row["dual_graph"]
        assert dual["component_count"] == 1
        assert dual["degree_min"] == dual["degree_max"] == 4
        assert dual["all_faces_two_sided_and_colour_consistent"] is True
        assert len(dual["graph_sha256"]) == 64
        assert all(edge["route_label"] in ("00", "01", "10", "11") for edge in dual["edges"])


def test_q2_pl_dual_refinement_records_parent_tetrahedra():
    result = run(max_generation=2)
    for row in result["generations"][1:]:
        assert all(
            len(vertex["flag_chain"]) == 4
            and len(vertex["parent_tetrahedron"]) == 4
            for vertex in row["dual_graph"]["vertices"]
        )
