import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(HERE))

import bqg_w5_sparse_active_path_support as SAP
import bqg_w5_su2_structural_basis as SB
import bqg_w5_pairwise_factor_graph as PG


def tiny_seed():
    spins = (0,) * 10
    labels = (("C", 0, 0, 0, 0),) + (("G", 0),) * 4
    return spins, labels


def assert_exact_included(word, source=0, jmax=3):
    key = tiny_seed()
    exact = SB.structural_sequence_support(
        key, word, source, jmax, charged_nodes=(0,), max_states=20000
    )
    sparse = SAP.sequence(
        PG._initial(key), word, source, jmax, charged_nodes=(0,), limit=200000
    )
    assert exact
    assert all(SAP.contains_joint(sparse, state) for state in exact)
    # The active relation count is strictly below the ten-pair global model.
    assert all(len(state["pairs"]) < 10 for state in sparse.values())


def test_sparse_active_path_includes_composite_hit_volume_hit_support():
    word = (("P", (0, 1, 2, 1, 0)), ("V", 0), ("P", (0, 2, 3, 2, 0)))
    assert_exact_included(word)


def test_sparse_active_path_includes_volume_inserted_before_hit():
    word = (("V", 0), ("P", (0, 1, 2, 0)))
    assert_exact_included(word)


def test_active_scope_is_local_to_primitive_paths():
    word = (("P", (0, 1, 2, 0)), ("V", 0))
    assert SAP.active_scope(word, 0) == (0, 1, 2)
    seed = PG._initial(tiny_seed())
    result = SAP.sequence(seed, word, 0, 3, charged_nodes=SAP.NODES,
                            limit=200000)
    prof = SAP.profile(result)
    assert prof["joint_keys_enumerated"] is False
    assert prof["amplitudes_evaluated"] is False
    assert prof["missing_pairs_are_unconstrained_cartesian_products"] is True
    terminal = SAP.project_final([result], (0,), limit=200000)
    assert terminal
    assert all(len(state["pairs"]) < 10 for branch in terminal
               for state in branch.values())


def test_sparse_sequential_words_keep_intermediate_nodes_charged():
    key = tiny_seed()
    first = (("P", (0, 1, 2, 1, 0)),)
    second = (("V", 0), ("P", (0, 2, 3, 2, 0)))
    exact = SB.structural_sequence_support(
        key, first + second, 0, 3, charged_nodes=(0,), max_states=20000
    )
    intermediate = SAP.sequence(
        PG._initial(key), first, 0, 3, charged_nodes=SAP.NODES,
        limit=200000
    )
    sparse = SAP.sequence(
        intermediate, second, 0, 3, charged_nodes=SAP.NODES,
        limit=200000
    )
    final = SAP.project_final([sparse], (0,), limit=200000)[0]
    assert exact
    assert all(SAP.contains_joint(final, state) for state in exact)


def test_he_spec_branches_preserve_inherited_nontrivial_pair_scope(monkeypatch):
    import k5_peter_weyl_safe_hda_column as PW

    left = ("C", 0, 0, 0, 0)
    left_alt = ("C", 2, 0, 0, 0)
    unary = tuple(
        frozenset((left, left_alt)) if v in (3, 4)
        else frozenset((left,))
        for v in range(5)
    )
    # A genuine inherited correlation on nodes outside this primitive word.
    relation = frozenset(((left, left),))
    incoming = {(0,) * 10: {
        "unary": unary,
        "pairs": {(3, 4): relation},
    }}

    monkeypatch.setattr(PW, "oriented_specs",
                        lambda source: ((1, (0, 1, 2, 3)),))
    monkeypatch.setattr(PW, "T_sequences",
                        lambda *args: ((1, (("V", 0),)),))
    monkeypatch.setattr(PW, "adjoint_sequence", lambda word: word)
    monkeypatch.setattr(SAP, "sequence",
                        lambda support, *args, **kwargs: support)

    result = SAP.he_spec_branches([incoming], 0, 3, charged_nodes=SAP.NODES)
    assert len(result) == 1
    state = result[0][(0,) * 10]
    assert (3, 4) in state["pairs"]
    assert state["pairs"][(3, 4)] == relation


def test_sparse_sequence_never_materializes_omitted_pair_factors(monkeypatch):
    def forbidden_expand(_support):
        raise AssertionError("sparse propagation must not expand all ten pairs")

    monkeypatch.setattr(SAP, "expand", forbidden_expand)
    result = SAP.sequence(
        PG._initial(tiny_seed()),
        (("P", (0, 1, 2, 0)), ("V", 0)),
        0, 3, charged_nodes=SAP.NODES, limit=200000
    )
    assert result
    assert all(len(state["pairs"]) < 10 for state in result.values())
