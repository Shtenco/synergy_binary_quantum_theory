"""Explicit Stinespring dilation of a two-outcome source measurement.

A source qubit is measured in an arbitrary orthonormal basis. A meter qubit
starts in |0>; a basis-controlled NOT correlates the source eigenstate with
the meter outcome. Reading the meter in its computational basis realizes
the projective instrument. All pre-existing environment registers remain
explicit until meter readout.

This is measurement theory, not a gravity model.
"""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

ComplexArray = NDArray[np.complex128]
__all__ = ["apply_source_measurement_dilation", "conditional_measurement_outcome", "measurement_outcome_probabilities", "source_measurement_dilation_unitary"]

def _validated_basis(basis: np.ndarray) -> ComplexArray:
    value = np.asarray(basis, dtype=np.complex128)
    if value.shape != (2, 2) or not np.allclose(value.conj().T @ value, np.eye(2), atol=1e-10):
        raise ValueError("measurement_basis must be a 2x2 unitary")
    return value

def source_measurement_dilation_unitary(measurement_basis: np.ndarray) -> ComplexArray:
    """Return U_SM = sum_k |b_k><b_k| tensor X^k (source, meter order)."""
    basis = _validated_basis(measurement_basis)
    identity = np.eye(2, dtype=np.complex128)
    bit_flip = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    p0 = np.outer(basis[:, 0], basis[:, 0].conj())
    p1 = np.outer(basis[:, 1], basis[:, 1].conj())
    return np.kron(p0, identity) + np.kron(p1, bit_flip)

def apply_source_measurement_dilation(state: np.ndarray, measurement_basis: np.ndarray) -> ComplexArray:
    """Append meter |0> and apply dilation; meter is the final tensor axis.

    Input shape is (2, *environment_shape). Output shape is
    (2, *environment_shape, 2). Meter slices are unnormalized outcome states.
    """
    value = np.asarray(state, dtype=np.complex128)
    if value.ndim < 1 or value.shape[0] != 2:
        raise ValueError("state must have source dimension 2 on axis zero")
    if not np.all(np.isfinite(value)):
        raise ValueError("state amplitudes must be finite")
    unitary = source_measurement_dilation_unitary(measurement_basis)
    env_shape = value.shape[1:]
    flat = value.reshape(2, -1)
    input_state = np.zeros((4, flat.shape[1]), dtype=np.complex128)
    input_state[0, :] = flat[0, :]
    input_state[2, :] = flat[1, :]
    evolved = (unitary @ input_state).reshape(2, 2, -1)
    return evolved.transpose(0, 2, 1).reshape((2,) + env_shape + (2,))

def measurement_outcome_probabilities(dilated_state: np.ndarray) -> tuple[float, float]:
    """Return probabilities encoded in the final two-level meter axis."""
    value = np.asarray(dilated_state, dtype=np.complex128)
    if value.ndim < 2 or value.shape[-1] != 2:
        raise ValueError("dilated state must have a final two-level meter axis")
    return tuple(float(np.sum(np.abs(value[..., k]) ** 2).real) for k in range(2))


def conditional_measurement_outcome(
    dilated_state: np.ndarray, outcome: int
) -> tuple[float, ComplexArray]:
    """Return (probability, normalized source+environment state) for one meter result.

    The meter axis is removed only after selecting its outcome. A zero-probability
    outcome raises ZeroDivisionError rather than returning an undefined state.
    """
    value = np.asarray(dilated_state, dtype=np.complex128)
    if value.ndim < 2 or value.shape[-1] != 2:
        raise ValueError("dilated state must have a final two-level meter axis")
    if outcome not in (0, 1):
        raise ValueError("outcome must be 0 or 1")
    branch = value[..., outcome]
    probability = float(np.sum(np.abs(branch) ** 2).real)
    if not np.isfinite(probability) or probability <= 1e-14:
        raise ZeroDivisionError("measurement outcome has zero or unresolved probability")
    return probability, branch / np.sqrt(probability)
