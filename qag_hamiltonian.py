"""Explicit Hamiltonian layer for the QAG toy-force benchmark.

The model is a one-dimensional probe in a branch-conditioned linear
potential V_R(x)=-F_R*x and V_L(x)=0. It is a controlled non-gravitational
toy model. On the finite periodic grid the linear potential is understood
only over the packet support; users must keep the packet far from the seam.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]


def linear_potential_hamiltonian(
    x: NDArray[np.float64], force: float
) -> ComplexMatrix:
    """Return H=V(X)=-F X, so the classical force is -dV/dx=+F."""
    grid = np.asarray(x, dtype=float)
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid)):
        raise ValueError("x must be a finite one-dimensional grid")
    if not np.isfinite(force):
        raise ValueError("force must be finite")
    return np.diag(-float(force) * grid).astype(np.complex128)


def unitary_from_diagonal_hamiltonian(
    hamiltonian: ComplexMatrix, duration: float, hbar: float = 1.0
) -> ComplexMatrix:
    """Compute exp(-i H t/hbar) for a validated diagonal Hermitian H."""
    h = np.asarray(hamiltonian, dtype=np.complex128)
    if h.ndim != 2 or h.shape[0] != h.shape[1]:
        raise ValueError("Hamiltonian must be square")
    if not np.allclose(h, h.conj().T, atol=1e-12):
        raise ValueError("Hamiltonian must be Hermitian")
    off_diagonal = h - np.diag(np.diag(h))
    if not np.allclose(off_diagonal, 0.0, atol=1e-12):
        raise ValueError("this exact helper accepts diagonal Hamiltonians only")
    if not np.isfinite(duration) or duration < 0.0:
        raise ValueError("duration must be finite and nonnegative")
    if not np.isfinite(hbar) or hbar <= 0.0:
        raise ValueError("hbar must be finite and positive")
    phases = np.exp(-1j * np.diag(h).real * duration / hbar)
    return np.diag(phases).astype(np.complex128)


def branch_unitaries_from_linear_potential(
    x: NDArray[np.float64],
    force_left: float,
    force_right: float,
    duration: float,
    hbar: float = 1.0,
) -> tuple[ComplexMatrix, ComplexMatrix]:
    """Build both controlled branch unitaries from explicit potential H's."""
    h_left = linear_potential_hamiltonian(x, force_left)
    h_right = linear_potential_hamiltonian(x, force_right)
    return (
        unitary_from_diagonal_hamiltonian(h_left, duration, hbar),
        unitary_from_diagonal_hamiltonian(h_right, duration, hbar),
    )
