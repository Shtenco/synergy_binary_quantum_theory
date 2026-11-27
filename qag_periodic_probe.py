"""Periodic-grid probe for exact momentum-kick postselection benchmarks.

This is a controlled toy interaction. A branch-dependent momentum kick is
not a Newtonian or relativistic gravitational Hamiltonian.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from qag_postselection import (
    conditional_state_via_kraus,
    expectation,
    weak_value_right,
)

ComplexMatrix = NDArray[np.complex128]


@dataclass(frozen=True)
class KickBenchmark:
    success_probability: float
    initial_momentum: float
    conditional_momentum: float
    conditional_shift: float
    weak_value_prediction: float
    relative_weak_error: float


def periodic_position_momentum_operators(
    n: int, length: float, hbar: float = 1.0
) -> tuple[NDArray[np.float64], ComplexMatrix]:
    """Return periodic position grid and Hermitian spectral momentum matrix."""
    if n < 8 or n % 2:
        raise ValueError("grid size n must be even and at least 8")
    if length <= 0.0 or hbar <= 0.0:
        raise ValueError("length and hbar must be positive")
    dx = length / n
    x = (np.arange(n, dtype=float) - n // 2) * dx
    wave_numbers = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    fourier = np.fft.fft(np.eye(n, dtype=np.complex128),
                         axis=0, norm="ortho")
    momentum = fourier.conj().T @ np.diag(hbar * wave_numbers) @ fourier
    momentum = (momentum + momentum.conj().T) / 2.0
    return x, momentum


def gaussian_probe_density(
    x: NDArray[np.float64], sigma: float, x0: float = 0.0,
    p0: float = 0.0, hbar: float = 1.0
) -> ComplexMatrix:
    """Normalized sampled Gaussian pure-state density matrix."""
    if sigma <= 0.0 or hbar <= 0.0:
        raise ValueError("sigma and hbar must be positive")
    amplitude = np.exp(-((x - x0) ** 2) / (4.0 * sigma**2))
    phase = np.exp(1j * p0 * x / hbar)
    state = np.asarray(amplitude * phase, dtype=np.complex128)
    norm = np.linalg.norm(state)
    if norm == 0.0:
        raise ValueError("sampled Gaussian has zero numerical norm")
    state /= norm
    return np.outer(state, state.conj())


def momentum_kick_unitary(
    x: NDArray[np.float64], delta_p: float, hbar: float = 1.0
) -> ComplexMatrix:
    """U=exp(i delta_p X/hbar), giving a positive delta_p momentum shift."""
    if hbar <= 0.0:
        raise ValueError("hbar must be positive")
    return np.diag(np.exp(1j * delta_p * x / hbar)).astype(np.complex128)


def run_weak_value_kick_benchmark(
    epsilon: float = 0.1,
    delta_p_right: float = -1.0e-3,
    n: int = 256,
    length: float = 40.0,
    sigma: float = 1.0,
    hbar: float = 1.0,
) -> KickBenchmark:
    """Compare exact postselected probe momentum with the weak-value limit.

    Branch L has no kick. Branch R receives delta_p_right. The source is
    prepared in the Bose–Vaidman two-branch parameterisation and postselected
    onto (|L>-|R>)/sqrt(2). The returned weak prediction is Re(A_w)*delta_p.
    """
    if not 0.0 < epsilon < 1.0:
        raise ValueError("epsilon must lie strictly between zero and one")
    x, momentum = periodic_position_momentum_operators(n, length, hbar)
    rho_probe = gaussian_probe_density(x, sigma, hbar=hbar)
    source = np.array(
        [np.sqrt(1.0 + epsilon), np.sqrt(1.0 - epsilon)],
        dtype=np.complex128,
    ) / np.sqrt(2.0)
    postselection = np.array([1.0, -1.0], dtype=np.complex128) / np.sqrt(2.0)
    u_left = np.eye(n, dtype=np.complex128)
    u_right = momentum_kick_unitary(x, delta_p_right, hbar)
    result = conditional_state_via_kraus(
        source, postselection, [u_left, u_right], rho_probe
    )
    p_initial = float(expectation(momentum, rho_probe).real)
    p_conditional = float(expectation(momentum, result.state).real)
    shift = p_conditional - p_initial
    weak_prediction = float(
        weak_value_right(source, postselection).real * delta_p_right
    )
    denominator = max(abs(weak_prediction), np.finfo(float).tiny)
    return KickBenchmark(
        success_probability=result.probability,
        initial_momentum=p_initial,
        conditional_momentum=p_conditional,
        conditional_shift=shift,
        weak_value_prediction=weak_prediction,
        relative_weak_error=abs(shift - weak_prediction) / denominator,
    )
