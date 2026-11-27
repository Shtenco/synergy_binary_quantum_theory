from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "research" / "depth6_recovered_source"))
sys.path.insert(0, str(ROOT / "research" / "depth6_no_prune"))

import bqg_block_apply_invariants as A  # noqa: E402


def test_initial_reduced_block_keeps_subthreshold_nonzero_rows_in_no_prune_mode():
    spins = (1,) * 10
    labels = A.gauss_basis_labels(spins)
    block = np.zeros((len(labels), 1), dtype=complex)
    block[0, 0] = 1.0e-13

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
    vector = np.array([1.0e-13], dtype=complex)
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
