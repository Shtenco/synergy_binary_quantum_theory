from scripts.q2_causal_pl_correspondence_gate import run


def test_q2_causal_to_pl_correspondence_preserves_adjacency_and_refinement():
    result = run(max_generation=2)
    assert result["passed"] is True
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
    assert "does not derive a physical measure" in result["claim_boundary"]



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
