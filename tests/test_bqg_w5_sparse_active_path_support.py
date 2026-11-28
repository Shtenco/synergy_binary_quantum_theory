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
    word = (("P", (0, 1, 2, 0)), ("V", 0), ("P", (0, 3, 4)))
    assert_exact_included(word)


def test_sparse_active_path_includes_volume_inserted_before_hit():
    word = (("V", 0), ("P", (0, 1, 2)))
    assert_exact_included(word)


def test_active_scope_is_local_to_primitive_paths():
    word = (("P", (0, 1, 2, 0)), ("V", 0))
    assert SAP.active_scope(word, 0) == (0, 1, 2)
    seed = PG._initial(tiny_seed())
    result = SAP.sequence(seed, word, 0, 3, limit=200000)
    prof = SAP.profile(result)
    assert prof["joint_keys_enumerated"] is False
    assert prof["amplitudes_evaluated"] is False
    assert prof["missing_pairs_are_unconstrained_cartesian_products"] is True
