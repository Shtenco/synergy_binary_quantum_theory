"""Unitary meter dilation for reading a finite-ring apparatus position.

The apparatus position label is copied into an N-level meter by modular
controlled addition. The pre-readout source/apparatus state is retained as an
explicit entangled state until the meter outcome is selected.

This is a finite-dimensional projective-readout model, not a calibrated
macroscopic detector.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

ComplexArray = NDArray[np.complex128]
__all__ = [
    "apparatus_readout_unitary",
    "apply_apparatus_readout_dilation",
    "conditional_apparatus_readout",
    "apparatus_readout_probabilities",
]


def _validate_state(state: np.ndarray) -> ComplexArray:
    value = np.asarray(state, dtype=np.complex128)
    if value.ndim != 2 or value.shape[0] != value.shape[1] or value.shape[0] < 2:
        raise ValueError("joint_state must be a square (N, N) source-apparatus array")
    if not np.all(np.isfinite(value)):
        raise ValueError("joint_state amplitudes must be finite")
    norm = float(np.vdot(value, value).real)
    if norm <= 0.0 or not np.isfinite(norm):
        raise ValueError("joint_state must have nonzero finite norm")
    return value / np.sqrt(norm)


def apparatus_readout_unitary(n: int) -> ComplexArray:
    """Return SUM |y,m> -> |y,m+y mod N> in apparatus-major order."""
    if not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError("n must be an integer >= 2")
    unitary = np.zeros((n * n, n * n), dtype=np.complex128)
    for y in range(n):
        for meter in range(n):
            source_index = y * n + meter
            target_index = y * n + ((meter + y) % n)
            unitary[target_index, source_index] = 1.0
    return unitary


def apply_apparatus_readout_dilation(joint_state: np.ndarray) -> ComplexArray:
    """Append |0> meter and coherently record apparatus position.

    Input shape: (N_source, N_apparatus). Output shape:
    (N_source, N_apparatus, N_meter). Meter outcome m selects apparatus
    position m, with the corresponding source-apparatus slice unnormalized.
    """
    value = _validate_state(joint_state)
    n = value.shape[0]
    dilated = np.zeros((n, n, n), dtype=np.complex128)
    for apparatus_position in range(n):
        dilated[:, apparatus_position, apparatus_position] = value[:, apparatus_position]
    return dilated


def apparatus_readout_probabilities(dilated_state: np.ndarray) -> tuple[float, ...]:
    """Return Born probabilities for all N meter outcomes."""
    value = np.asarray(dilated_state, dtype=np.complex128)
    if value.ndim != 3 or value.shape[0] != value.shape[1] or value.shape[1] != value.shape[2]:
        raise ValueError("dilated_state must have shape (N, N, N)")
    if not np.all(np.isfinite(value)):
        raise ValueError("dilated_state amplitudes must be finite")
    return tuple(float(np.sum(np.abs(value[:, :, k]) ** 2).real) for k in range(value.shape[2]))


def conditional_apparatus_readout(
    dilated_state: np.ndarray, outcome: int
) -> tuple[float, ComplexArray]:
    """Return (probability, normalized source-apparatus state) for meter outcome."""
    value = np.asarray(dilated_state, dtype=np.complex128)
    if value.ndim != 3 or value.shape[0] != value.shape[1] or value.shape[1] != value.shape[2]:
        raise ValueError("dilated_state must have shape (N, N, N)")
    n = value.shape[2]
    if not isinstance(outcome, (int, np.integer)) or not 0 <= outcome < n:
        raise ValueError("outcome must be an integer in [0, N)")
    branch = value[:, :, int(outcome)]
    probability = float(np.vdot(branch, branch).real)
    if not np.isfinite(probability) or probability <= 1e-14:
        raise ZeroDivisionError("apparatus readout outcome has zero or unresolved probability")
    return probability, branch / np.sqrt(probability)
