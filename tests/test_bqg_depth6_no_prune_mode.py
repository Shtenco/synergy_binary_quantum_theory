from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "research" / "depth6_recovered_source"))
sys.path.insert(0, str(ROOT / "research" / "depth6_no_prune"))

import bqg_block_apply_invariants as A  # noqa: E402
import bqg_depth6_generic_master_row as M  # noqa: E402
import github_mixed_master_shard as Runner  # noqa: E402


def test_initial_reduced_block_keeps_subthreshold_nonzero_rows_in_no_prune_mode():
    spins = (1,) * 10
    labels = A.gauss_basis_labels(spins)
    block = np.zeros((len(labels), 1), dtype=complex)
    block[0, 0] = 1.0e-200

    legacy = A.initial_reduced_block(spins, block)
    no_prune = A.initial_reduced_block(spins, block, tol=0.0)

    assert legacy == {}
    assert len(no_prune) == 1
    np.testing.assert_array_equal(next(iter(no_prune.values())), block[0])


def test_sparse_accumulation_is_order_independent_when_pruning_is_disabled():
    tau = 1.0e-12
    forward = (1.2 * tau, -0.7 * tau, 0.7 * tau)
    reverse = tuple(reversed(forward))

    def accumulate(values, tol):
        result = {}
        for value in values:
            A.add_vec(result, "same-key", value, tol=tol)
        return result.get("same-key", 0.0)

    assert accumulate(forward, tol=0.0) == accumulate(reverse, tol=0.0)
    assert np.isclose(accumulate(forward, tol=0.0), 1.2 * tau, rtol=0.0, atol=1e-27)
    assert accumulate(forward, tol=tau) != accumulate(reverse, tol=tau)


def test_gauss_projection_zero_tolerance_retains_subthreshold_matrix_entries():
    spins = (1,) * 10
    basis_label = A.gauss_basis_labels(spins)[0]
    locs = tuple((K, K, 0, 0) for K in basis_label)
    key = (spins, locs)
    vector = np.array([1.0e-200], dtype=complex)
    reduced = {key: vector}

    legacy = A.project_gauss_block(reduced)
    no_prune = A.project_gauss_block(reduced, tol=0.0)

    assert legacy == {}
    assert spins in no_prune
    assert no_prune[spins][0, 0] == vector[0]


def test_zero_outer_tolerance_propagates_through_the_complete_node_wrapper(monkeypatch):
    observed = {}

    def initial(spins, block, tol=A.TOL):
        observed["initial_tol"] = tol
        return {}

    def oriented(state, vertex, a, b, c, jmax, tol):
        observed["oriented_tol"] = tol
        return {}

    def project(state, tol):
        observed["project_tol"] = tol
        return {}

    monkeypatch.setattr(A, "initial_reduced_block", initial)
    monkeypatch.setattr(A.R, "oriented_specs", lambda vertex: ((1, (vertex, 1, 2, 3)),))
    monkeypatch.setattr(A, "oriented_T_ultrafast", oriented)
    monkeypatch.setattr(A, "project_gauss_block", project)

    A.H_node_apply_ultrafast((1,) * 10, 0, np.zeros((32, 1)), tol=0.0)

    assert observed == {
        "initial_tol": 0.0,
        "oriented_tol": 0.0,
        "project_tol": 0.0,
    }



def test_no_prune_master_row_resolves_only_isolated_derived_operator_modules():
    no_prune_root = ROOT / "research" / "depth6_no_prune"
    assert Path(M.A.__file__).resolve().parent == no_prune_root
    assert Path(M.A.H.__file__).resolve().parent == no_prune_root
    assert Path(M.A.Rec.__file__).resolve().parent == no_prune_root
    assert Path(M.J.__file__).resolve().parent == no_prune_root
    assert M.J.expand.__defaults__[-1] == 0.0
    assert M.J.fast_selector.__defaults__[-1] == 1e-10
    assert M.master_map.__defaults__[-1] == 0.0
    assert M.reduced_s3_node.__defaults__[-1] == 0.0
    assert M.raw_to_s2_coords.__defaults__[-1] == 0.0
    assert M.reduced_s2_node.__defaults__[-1] == 0.0



def test_no_prune_shard_runner_uses_the_isolated_master_and_explicit_zero_tolerance(
    monkeypatch, tmp_path
):
    observed = {}
    assert Runner.M is M
    assert Runner.J is M.J

    monkeypatch.setattr(Runner, "_RECS", [(7, (1,) * 10, 1, 1, 1.0)])
    monkeypatch.setattr(Runner, "_KEY", "32")
    monkeypatch.setattr(Runner, "_OUT", tmp_path)
    monkeypatch.setattr(Runner.J, "fast_selector", lambda rep, key: np.ones((1, 1)))

    def master_map(rep, key, selector, tol):
        observed["tol"] = tol
        return {}, selector

    monkeypatch.setattr(Runner.M, "master_map", master_map)
    result = Runner.work(0)

    assert result["ok"] is True
    assert observed["tol"] == 0.0
    assert (tmp_path / "b00007.pkl").is_file()



def test_no_prune_clebsch_gordan_reconstruction_keeps_subthreshold_coefficients(
    monkeypatch,
):
    recoupling = A.Rec
    recoupling.coupled_tensor.cache_clear()
    monkeypatch.setattr(recoupling.R, "cg2", lambda *args: 1.0e-16)

    tensor = recoupling.coupled_tensor((1, 1, 1, 1), 0, 0, 0, 0)

    assert np.count_nonzero(tensor) > 0
    assert np.max(np.abs(tensor)) < 1.0e-40
    recoupling.coupled_tensor.cache_clear()



def test_no_prune_selector_expansion_retains_subthreshold_coordinate_segments(monkeypatch):
    selector = M.J
    monkeypatch.setattr(
        selector,
        "layout",
        lambda rep, seed: (("q", "zero"), {"q": slice(0, 1), "zero": slice(1, 2)}, 2),
    )
    monkeypatch.setattr(selector, "seed_basis", lambda q, seed: np.ones((1, 1)))
    monkeypatch.setattr(selector, "S3", ((0, 1, 2, 3, 4),))
    monkeypatch.setattr(selector.Z, "s3_stabilizer", lambda q: tuple(range(6)))
    monkeypatch.setattr(
        selector.X,
        "transform_coeff",
        lambda q, permutation, block: (q, block),
    )

    tiny = 1.0e-14
    expanded = selector.expand(
        (1,) * 10, "triv", np.array([[tiny], [0.0]], dtype=complex)
    )

    assert tuple(expanded) == ("q",)
    assert expanded["q"][0, 0] == tiny



def test_zero_tolerance_retains_values_whose_squared_norm_underflows():
    tiny = np.array([1.0e-200 + 0.0j])
    assert np.linalg.norm(tiny) == 0.0
    assert A._keep_nonzero_or_above(tiny, 0.0)
    assert M._keep_nonzero_or_above(tiny, 0.0)
    assert not A._keep_nonzero_or_above(tiny, 1.0e-12)



def test_small_real_nonvacuum_master_block_runs_end_to_end_without_pruning():
    # A recovered depth-6 orbit representative with a small selector coordinate
    # dimension. This executes the actual selector and master map, not a mock.
    rep = (1, 1, 1, 1, 1, 1, 3, 3, 1, 5)
    key = "32"
    expected_multiplicity = M.J.exact_mult(rep, key)
    assert expected_multiplicity > 0

    W = M.J.fast_selector(rep, key)
    assert W.shape == (M.J.layout(rep, M.J.CFG[key][0])[2], expected_multiplicity)

    master_rows, returned_W = M.master_map(rep, key, W, tol=0.0)

    assert returned_W.shape == W.shape
    assert master_rows
    for matrix in master_rows.values():
        assert matrix.ndim == 2
        assert matrix.shape[1] == expected_multiplicity
        assert np.isfinite(matrix).all()
