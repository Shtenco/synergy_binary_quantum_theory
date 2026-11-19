import numpy as np
import pytest
from scipy.sparse import csr_matrix

from scripts.q2_dirichlet_form_convergence import (
    GraphDirichletForm,
    recovery_probe_report,
    validate_prolongation,
)


def _uniform_interval_form(cells):
    h = 1.0 / cells
    edges = np.column_stack((np.arange(cells), np.arange(1, cells + 1)))
    weights = np.full(cells, 1.0 / h)
    measure = np.full(cells + 1, h)
    measure[[0, -1]] *= 0.5
    return GraphDirichletForm.build(cells + 1, edges, weights, measure)


def test_interval_linear_function_has_exact_dirichlet_energy():
    form = _uniform_interval_form(32)
    x = np.linspace(0.0, 1.0, 33)
    assert form.energy(x) == pytest.approx(1.0, abs=1e-13)
    assert form.mass_norm_squared(x) == pytest.approx(1.0 / 3.0, abs=2e-4)


def test_graph_form_generator_annihilates_constants_and_is_symmetric():
    form = _uniform_interval_form(12)
    assert np.allclose(form.apply_generator(np.ones(13)), 0.0, atol=1e-13)
    assert np.allclose(form.stiffness.toarray(), form.stiffness.toarray().T)
    assert np.linalg.eigvalsh(form.stiffness.toarray())[0] >= -1e-12
    assert form.summary()["measure_total"] == pytest.approx(1.0)


def test_mass_weighted_resolvent_has_small_residual():
    form = _uniform_interval_form(20)
    rhs = np.sin(np.linspace(0.0, np.pi, 21))
    solution, residual = form.resolvent_action(rhs, z=0.7)
    assert np.all(np.isfinite(solution))
    assert residual < 1e-12


def test_prolongation_must_preserve_constants_and_be_nonnegative():
    p = csr_matrix([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]])
    checked, error = validate_prolongation(p, n_fine=3, n_coarse=2)
    assert checked.shape == (3, 2)
    assert error == 0.0
    with pytest.raises(ValueError, match="preserve constants"):
        validate_prolongation(csr_matrix([[0.5, 0.0]]), n_fine=1, n_coarse=2)
    with pytest.raises(ValueError, match="nonnegative"):
        validate_prolongation(csr_matrix([[1.2, -0.2]]), n_fine=1, n_coarse=2)


def test_recovery_probe_report_is_explicitly_not_a_mosco_proof():
    coarse = GraphDirichletForm.build(
        2, np.array([[0, 1]]), np.array([1.0]), np.array([0.5, 0.5])
    )
    fine = GraphDirichletForm.build(
        3, np.array([[0, 1], [1, 2]]), np.array([2.0, 2.0]),
        np.array([0.25, 0.5, 0.25])
    )
    prolongation = csr_matrix([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]])
    report = recovery_probe_report(
        coarse, fine, prolongation, [np.array([0.0, 1.0]), np.ones(2)]
    )
    assert report["status"] == "FINITE_RECOVERY_PROBES_ONLY"
    assert len(report["probes"]) == 2
    assert "no liminf proof" in report["scope_boundary"]
