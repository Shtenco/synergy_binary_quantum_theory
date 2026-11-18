from scripts.q2_causal_graph_export import build, expected_vertices


def test_q2_causal_graph_matches_exact_vertex_and_edge_growth():
    result = build(max_generation=5)
    assert result["passed"] is True
    assert [
        (row["vertex_count"], row["active_edge_count"])
        for row in result["generations"]
    ] == [
        (2, 1),
        (6, 8),
        (38, 64),
        (294, 512),
        (2342, 4096),
        (18726, 32768),
    ]
    assert all(
        row["vertex_count"] == expected_vertices(row["generation"])
        for row in result["generations"]
    )


def test_q2_causal_graph_is_connected_dag_with_exact_causal_depth():
    result = build(max_generation=5)
    for row in result["generations"]:
        assert row["connected"] is True
        assert row["directed_acyclic"] is True
        assert row["path_length"] == 2 ** row["generation"]
        assert row["source_to_sink_path_count"] == 4 ** (2 ** row["generation"] - 1)


def test_q2_causal_graph_preserves_route_vertex_ancestry():
    result = build(max_generation=3)
    for generation, row in enumerate(result["generations"][:-1], start=0):
        mapping = row["refinement_map_to_next"]
        assert len(mapping) == 8 ** generation
        for parent in mapping:
            assert len(parent["route_vertices"]) == 4
            assert len(parent["child_edges"]) == 8
            assert [x["label"] for x in parent["route_vertices"]] == [
                "00", "01", "10", "11"
            ]
    assert result["pl_identification"] is None
    assert result["qfi_state_assignment"] is None
