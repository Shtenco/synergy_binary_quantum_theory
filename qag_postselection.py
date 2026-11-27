"""Finite-dimensional postselection engine for QAG-1.

This module models a controlled interaction between a two-branch source and
a finite-dimensional probe. It is a quantum-information benchmark, not a
gravity simulator and not evidence for physical antigravity.

Dependencies: NumPy.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]
ComplexVector = NDArray[np.complex128]


def _as_complex_matrix(value: np.ndarray, name: str) -> ComplexMatrix:
    matrix = np.asarray(value, dtype=np.complex128)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    return matrix


def _validate_density(rho: np.ndarray, name: str = "density matrix",
                      atol: float = 1e-10) -> ComplexMatrix:
    matrix = _as_complex_matrix(rho, name)
    if not np.allclose(matrix, matrix.conj().T, atol=atol):
        raise ValueError(f"{name} must be Hermitian")
    if not np.isclose(np.trace(matrix), 1.0, atol=atol):
        raise ValueError(f"{name} must have unit trace")
    eig = np.linalg.eigvalsh(matrix)
    if eig.min(initial=0.0) < -atol:
        raise ValueError(f"{name} must be positive semidefinite")
    return matrix


def _validate_unitary(unitary: np.ndarray, dimension: int,
                      name: str, atol: float = 1e-10) -> ComplexMatrix:
    matrix = _as_complex_matrix(unitary, name)
    if matrix.shape != (dimension, dimension):
        raise ValueError(f"{name} has incompatible dimension")
    if not np.allclose(matrix.conj().T @ matrix,
                       np.eye(dimension), atol=atol):
        raise ValueError(f"{name} must be unitary")
    return matrix


@dataclass(frozen=True)
class ConditionalResult:
    """Unnormalised and normalised postselected probe state."""

    probability: float
    unnormalised: ComplexMatrix
    state: ComplexMatrix


def source_density_from_amplitudes(amplitudes: Sequence[complex],
                                   coherence: complex = 1.0) -> ComplexMatrix:
    """Two-branch source density matrix with controllable off-diagonal coherence.

    coherence=1 gives the pure state |a><a|; coherence=0 gives its
    fully dephased branch mixture. Positivity requires |coherence| <= 1.
    """
    a = np.asarray(amplitudes, dtype=np.complex128)
    if a.shape != (2,):
        raise ValueError("exactly two source amplitudes are required")
    norm = float(np.vdot(a, a).real)
    if norm <= 0.0:
        raise ValueError("source amplitudes must have nonzero norm")
    a = a / np.sqrt(norm)
    nu = complex(coherence)
    if abs(nu) > 1.0 + 1e-12:
        raise ValueError("absolute coherence must not exceed one")
    rho = np.outer(a, a.conj())
    rho[0, 1] *= nu
    rho[1, 0] *= nu.conjugate()
    return _validate_density(rho, "source density")


def postselection_kraus(source_amplitudes: Sequence[complex],
                        postselection: Sequence[complex],
                        branch_unitaries: Sequence[np.ndarray]) -> ComplexMatrix:
    """Return K=<phi|U|psi> for a pure source and rank-one postselection."""
    a = np.asarray(source_amplitudes, dtype=np.complex128)
    phi = np.asarray(postselection, dtype=np.complex128)
    if a.shape != (2,) or phi.shape != (2,):
        raise ValueError("source and postselection must be two-component vectors")
    if not np.isclose(np.vdot(a, a).real, 1.0, atol=1e-10):
        raise ValueError("source amplitudes must be normalised")
    if not np.isclose(np.vdot(phi, phi).real, 1.0, atol=1e-10):
        raise ValueError("postselection vector must be normalised")
    if len(branch_unitaries) != 2:
        raise ValueError("exactly two branch unitaries are required")
    u0 = _as_complex_matrix(branch_unitaries[0], "U_L")
    u1 = _validate_unitary(branch_unitaries[1], u0.shape[0], "U_R")
    u0 = _validate_unitary(u0, u0.shape[0], "U_L")
    return np.conjugate(phi[0]) * a[0] * u0 + np.conjugate(phi[1]) * a[1] * u1


def conditional_probe_state(source_density: np.ndarray,
                             postselection: Sequence[complex],
                             branch_unitaries: Sequence[np.ndarray],
                             probe_density: np.ndarray,
                             atol: float = 1e-12) -> ConditionalResult:
    """Evaluate the exact conditional channel for a possibly mixed source.

    rho_tilde = sum_ij rho_source[i,j] conj(phi[i]) phi[j]
                         U_i rho_probe U_j^dagger.

    The source and postselection are two-dimensional; the probe dimension is
    arbitrary. The result is not renormalised until its success probability
    has been computed.
    """
    rho_s = _validate_density(source_density, "source density")
    rho_p = _validate_density(probe_density, "probe density")
    phi = np.asarray(postselection, dtype=np.complex128)
    if phi.shape != (2,) or not np.isclose(np.vdot(phi, phi).real, 1.0,
                                           atol=1e-10):
        raise ValueError("postselection must be a normalised two-vector")
    if len(branch_unitaries) != 2:
        raise ValueError("exactly two branch unitaries are required")
    dim = rho_p.shape[0]
    unitaries = [
        _validate_unitary(branch_unitaries[i], dim, f"U_{i}")
        for i in range(2)
    ]
    out = np.zeros_like(rho_p)
    for i in range(2):
        for j in range(2):
            coefficient = rho_s[i, j] * np.conjugate(phi[i]) * phi[j]
            out += coefficient * unitaries[i] @ rho_p @ unitaries[j].conj().T
    out = (out + out.conj().T) / 2.0
    probability_complex = np.trace(out)
    if abs(probability_complex.imag) > 1e-9:
        raise ArithmeticError("postselection probability acquired an imaginary part")
    probability = float(probability_complex.real)
    if probability < -atol or probability > 1.0 + 1e-9:
        raise ArithmeticError(f"invalid postselection probability: {probability}")
    if probability <= atol:
        raise ZeroDivisionError("postselection probability is zero or numerically unresolved")
    state = out / probability
    _validate_density(state, "conditional probe state", atol=1e-8)
    return ConditionalResult(probability, out, state)


def conditional_state_via_kraus(source_amplitudes: Sequence[complex],
                                postselection: Sequence[complex],
                                branch_unitaries: Sequence[np.ndarray],
                                probe_density: np.ndarray) -> ConditionalResult:
    """Pure-source reference path using the single Kraus operator."""
    rho_p = _validate_density(probe_density, "probe density")
    k = postselection_kraus(source_amplitudes, postselection, branch_unitaries)
    unnormalised = k @ rho_p @ k.conj().T
    probability = float(np.trace(unnormalised).real)
    if probability <= 1e-12:
        raise ZeroDivisionError("postselection probability is zero or unresolved")
    state = unnormalised / probability
    return ConditionalResult(probability, unnormalised, state)


def complete_two_outcome_probabilities(source_density: np.ndarray,
                                       probe_density: np.ndarray,
                                       branch_unitaries: Sequence[np.ndarray],
                                       basis: np.ndarray | None = None) -> tuple[float, float]:
    """Return probabilities for an orthonormal two-outcome source measurement."""
    if basis is None:
        basis = np.eye(2, dtype=np.complex128)
    basis = np.asarray(basis, dtype=np.complex128)
    if basis.shape != (2, 2) or not np.allclose(
            basis.conj().T @ basis, np.eye(2), atol=1e-10):
        raise ValueError("measurement basis must be a 2x2 unitary matrix")
    rho_s = _validate_density(source_density, "source density")
    rho_p = _validate_density(probe_density, "probe density")
    dim = rho_p.shape[0]
    unitaries = [
        _validate_unitary(branch_unitaries[i], dim, f"U_{i}")
        for i in range(2)
    ]
    probs = []
    for column in range(2):
        phi = basis[:, column]
        unnormalised = np.zeros_like(rho_p)
        for i in range(2):
            for j in range(2):
                unnormalised += (
                    rho_s[i, j] * np.conjugate(phi[i]) * phi[j]
                    * unitaries[i] @ rho_p @ unitaries[j].conj().T
                )
        probs.append(float(np.trace(unnormalised).real))
    if not np.isclose(sum(probs), 1.0, atol=1e-9):
        raise ArithmeticError("complete measurement probabilities do not sum to one")
    return float(probs[0]), float(probs[1])


def weak_value_right(source_amplitudes: Sequence[complex],
                     postselection: Sequence[complex]) -> complex:
    """Weak value of |R><R|; raises at exactly orthogonal pre/postselection."""
    a = np.asarray(source_amplitudes, dtype=np.complex128)
    phi = np.asarray(postselection, dtype=np.complex128)
    if a.shape != (2,) or phi.shape != (2,):
        raise ValueError("two-component vectors required")
    denominator = np.vdot(phi, a)
    if abs(denominator) < 1e-14:
        raise ZeroDivisionError("pre- and postselected states are orthogonal")
    return np.conjugate(phi[1]) * a[1] / denominator


def pauli_y_rotation(theta: float) -> ComplexMatrix:
    """Unitary exp(-i theta sigma_y/2), useful for a controlled toy benchmark."""
    c, s = np.cos(theta / 2.0), np.sin(theta / 2.0)
    return np.array([[c, -s], [s, c]], dtype=np.complex128)


def expectation(operator: np.ndarray, density: np.ndarray) -> complex:
    """Expectation value Tr(operator rho), with shape validation."""
    op = _as_complex_matrix(operator, "operator")
    rho = _validate_density(density)
    if op.shape != rho.shape:
        raise ValueError("operator and density dimensions differ")
    return complex(np.trace(op @ rho))
