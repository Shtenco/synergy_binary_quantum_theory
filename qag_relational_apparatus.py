"""Relational source-apparatus interaction on a finite cyclic position lattice.

The interaction is diagonal in the relative coordinate x_source-x_apparatus.
It therefore commutes exactly with simultaneous cyclic translations. The
conserved symmetry is a finite-ring translation operator (equivalently a
modular momentum sector), not an unbounded canonical momentum.

This is a quantum-information toy model. It is not a calibrated detector,
continuous-space mechanics, or a gravitational interaction.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray

ComplexArray = NDArray[np.complex128]
__all__ = [
    "RelationalInteractionResult",
    "apply_relational_interaction",
    "cyclic_translation",
    "gaussian_position_state",
    "relative_interaction_unitary",
    "translation_commutator_norm",
]


def _normalise(state: np.ndarray, n: int, name: str) -> ComplexArray:
    value = np.asarray(state, dtype=np.complex128)
    if value.shape != (n,):
        raise ValueError(f"{name} must have shape ({n},)")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must contain finite amplitudes")
    norm = float(np.vdot(value, value).real)
    if norm <= 0.0 or not np.isfinite(norm):
        raise ValueError(f"{name} must have nonzero finite norm")
    return value / np.sqrt(norm)


def gaussian_position_state(
    n: int, center: float, width: float
) -> ComplexArray:
    """Periodic Gaussian packet sampled on sites 0,...,n-1."""
    if not isinstance(n, (int, np.integer)) or n < 4:
        raise ValueError("n must be an integer >= 4")
    if not np.isfinite(center) or not np.isfinite(width) or width <= 0.0:
        raise ValueError("center must be finite and width positive")
    sites = np.arange(n, dtype=float)
    delta = (sites - float(center) + n / 2.0) % n - n / 2.0
    values = np.exp(-0.5 * (delta / width) ** 2).astype(np.complex128)
    return _normalise(values, n, "Gaussian position state")


def cyclic_translation(n: int, steps: int = 1) -> ComplexArray:
    """Matrix mapping |x> to |x+steps mod n>."""
    if not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError("n must be an integer >= 2")
    if not isinstance(steps, (int, np.integer)):
        raise ValueError("steps must be an integer")
    return np.roll(np.eye(n, dtype=np.complex128), int(steps) % n, axis=0)


def relative_interaction_unitary(
    n: int,
    coupling: float,
    duration: float,
    hbar: float = 1.0,
    harmonic: int = 1,
) -> ComplexArray:
    """Return U=exp[-i tau V(x_s-x_a)/hbar] in source-major basis.

    V(r)=coupling*cos(2*pi*harmonic*r/n). The periodic potential is a
    deliberately simple relational coupling; its sign and scale are model
    parameters, not gravitational constants.
    """
    if not isinstance(n, (int, np.integer)) or n < 4:
        raise ValueError("n must be an integer >= 4")
    if not isinstance(harmonic, (int, np.integer)) or harmonic < 1:
        raise ValueError("harmonic must be a positive integer")
    if not all(np.isfinite(v) for v in (coupling, duration, hbar)) or hbar <= 0:
        raise ValueError("coupling and duration must be finite; hbar positive")
    relative = (np.arange(n)[:, None] - np.arange(n)[None, :]) % n
    potential = coupling * np.cos(2.0 * np.pi * harmonic * relative / n)
    phases = np.exp(-1j * duration * potential / hbar).reshape(n * n)
    return np.diag(phases.astype(np.complex128))


def _entropy(rho: ComplexArray) -> float:
    hermitian = (rho + rho.conj().T) / 2.0
    eigenvalues = np.linalg.eigvalsh(hermitian).real
    eigenvalues = np.clip(eigenvalues, 0.0, 1.0)
    positive = eigenvalues[eigenvalues > 1e-14]
    return float(-np.sum(positive * np.log2(positive)))


@dataclass(frozen=True)
class RelationalInteractionResult:
    joint_state: ComplexArray
    joint_position_probability: NDArray[np.float64]
    source_density: ComplexArray
    apparatus_density: ComplexArray
    entanglement_entropy_bits: float
    translation_expectation_initial: complex
    translation_expectation_final: complex
    commutator_norm: float
    norm: float


def translation_commutator_norm(unitary: np.ndarray, n: int) -> float:
    """Frobenius norm of [U,T_source tensor T_apparatus]."""
    value = np.asarray(unitary, dtype=np.complex128)
    if value.shape != (n * n, n * n):
        raise ValueError("unitary must have shape (n*n, n*n)")
    translation = np.kron(cyclic_translation(n), cyclic_translation(n))
    return float(np.linalg.norm(value @ translation - translation @ value))


def apply_relational_interaction(
    source: np.ndarray,
    apparatus: np.ndarray,
    coupling: float,
    duration: float,
    hbar: float = 1.0,
    harmonic: int = 1,
) -> RelationalInteractionResult:
    """Entangle source and apparatus through a translation-invariant phase.

    Tensor order is source position first, apparatus position second. The
    returned joint position probabilities are also the direct readout table.
    """
    source_value = np.asarray(source, dtype=np.complex128)
    apparatus_value = np.asarray(apparatus, dtype=np.complex128)
    if source_value.ndim != 1 or apparatus_value.ndim != 1:
        raise ValueError("source and apparatus must be one-dimensional states")
    if source_value.size != apparatus_value.size:
        raise ValueError("source and apparatus dimensions must match")
    n = source_value.size
    source_value = _normalise(source_value, n, "source")
    apparatus_value = _normalise(apparatus_value, n, "apparatus")
    initial = np.outer(source_value, apparatus_value)
    unitary = relative_interaction_unitary(n, coupling, duration, hbar, harmonic)
    final = (unitary @ initial.reshape(n * n)).reshape(n, n)

    source_density = np.einsum("xy,zy->xz", final, final.conj())
    apparatus_density = np.einsum("xy,xz->yz", final, final.conj())
    source_density = (source_density + source_density.conj().T) / 2.0
    apparatus_density = (apparatus_density + apparatus_density.conj().T) / 2.0

    translation = np.kron(cyclic_translation(n), cyclic_translation(n))
    initial_translation = complex(np.vdot(initial.reshape(-1), translation @ initial.reshape(-1)))
    final_translation = complex(np.vdot(final.reshape(-1), translation @ final.reshape(-1)))
    norm = float(np.vdot(final.reshape(-1), final.reshape(-1)).real)
    entropy = _entropy(source_density)
    return RelationalInteractionResult(
        joint_state=final,
        joint_position_probability=np.abs(final) ** 2,
        source_density=source_density,
        apparatus_density=apparatus_density,
        entanglement_entropy_bits=entropy,
        translation_expectation_initial=initial_translation,
        translation_expectation_final=final_translation,
        commutator_norm=translation_commutator_norm(unitary, n),
        norm=norm,
    )
