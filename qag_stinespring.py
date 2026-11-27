"""Four-register Stinespring toy model with recoil and a measurement pointer.

Registers: source branch qubit, source recoil momentum lattice, probe momentum
lattice, and a neutral two-level pointer. Branch R transfers +q momentum units
to the probe and -q to the source recoil register. The pointer weakly records
the branch through a controlled rotation. Momentum conservation is exact
modulo the finite cyclic lattice; physical (non-wrapped) interpretation
requires negligible boundary probability.

This is a finite-dimensional quantum-information model, not gravity.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray

ComplexArray = NDArray[np.complex128]


def _normalised_vector(value: np.ndarray, dimension: int, name: str) -> ComplexArray:
    vector = np.asarray(value, dtype=np.complex128)
    if vector.shape != (dimension,):
        raise ValueError(f"{name} must have shape ({dimension},)")
    norm = float(np.vdot(vector, vector).real)
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError(f"{name} must have nonzero finite norm")
    return vector / np.sqrt(norm)


def gaussian_momentum_state(
    n: int, center: float = 0.0, width: float = 2.0
) -> ComplexArray:
    """Normalized Gaussian on centered integer momentum labels."""
    if n < 8 or n % 2 or width <= 0.0 or not np.isfinite(width):
        raise ValueError("n must be even >= 8 and width must be positive")
    labels = np.arange(n, dtype=float) - n // 2
    state = np.exp(-0.5 * ((labels - center) / width) ** 2).astype(complex)
    return _normalised_vector(state, n, "Gaussian momentum state")


def ry(theta: float) -> ComplexArray:
    """Real SU(2) rotation exp(-i theta sigma_y/2)."""
    if not np.isfinite(theta):
        raise ValueError("theta must be finite")
    c, s = np.cos(theta / 2.0), np.sin(theta / 2.0)
    return np.array([[c, -s], [s, c]], dtype=np.complex128)


@dataclass(frozen=True)
class StinespringResult:
    probability: float
    probe_state: ComplexArray
    joint_recoil_probe_pointer: ComplexArray
    total_momentum_initial: float
    total_momentum_conditional: float
    source_recoil_shift: float
    probe_shift: float
    pointer_branch_overlap: complex


def _momentum_expectation(state: ComplexArray, axis: int) -> float:
    n = state.shape[axis]
    labels = np.arange(n, dtype=float) - n // 2
    probabilities = np.sum(np.abs(state) ** 2, axis=tuple(
        i for i in range(state.ndim) if i != axis
    ))
    return float(np.dot(labels, probabilities).real)


def run_stinespring_recoil_model(
    source_amplitudes: np.ndarray,
    postselection: np.ndarray,
    source_recoil: np.ndarray,
    probe: np.ndarray,
    impulse_units: int = 1,
    pointer_angle: float = 0.0,
) -> StinespringResult:
    """Evolve all registers, postselect source branch, and trace ancillas.

    Branch L has no impulse and leaves the pointer unchanged. Branch R applies
    recoil/probe shifts (-q,+q) and rotates the pointer by pointer_angle.
    The initial pointer is |0>. The returned probe state is normalized after
    conditioning. All operations are unitary before postselection.
    """
    a = np.asarray(source_amplitudes, dtype=np.complex128)
    phi = np.asarray(postselection, dtype=np.complex128)
    if a.shape != (2,) or phi.shape != (2,):
        raise ValueError("source and postselection must be two-vectors")
    a = _normalised_vector(a, 2, "source amplitudes")
    phi = _normalised_vector(phi, 2, "postselection")
    recoil = np.asarray(source_recoil, dtype=np.complex128)
    probe_state = np.asarray(probe, dtype=np.complex128)
    if recoil.ndim != 1 or probe_state.ndim != 1 or recoil.size != probe_state.size:
        raise ValueError("source recoil and probe must be equal-size vectors")
    n = recoil.size
    if n < 8 or n % 2:
        raise ValueError("momentum lattice size must be even and >= 8")
    recoil = _normalised_vector(recoil, n, "source recoil")
    probe_state = _normalised_vector(probe_state, n, "probe")
    if not isinstance(impulse_units, (int, np.integer)):
        raise ValueError("impulse_units must be an integer")
    q = int(impulse_units)
    if abs(q) >= n // 2:
        raise ValueError("impulse magnitude must be smaller than half the lattice")

    pointer_initial = np.array([1.0, 0.0], dtype=np.complex128)
    pointer_r = ry(pointer_angle) @ pointer_initial
    branch_l = np.einsum("r,p,a->rpa", recoil, probe_state, pointer_initial)
    branch_r = np.einsum("r,p,a->rpa", np.roll(recoil, -q), 
                         np.roll(probe_state, q), pointer_r)
    joint = a[0] * branch_l + a[1] * branch_r

    initial_total = _momentum_expectation(recoil, 0) + _momentum_expectation(probe_state, 0)
    # Postselection contracts the source-branch index. Keep all physical ancillas.
    selected = np.conjugate(phi[0]) * a[0] * branch_l + np.conjugate(phi[1]) * a[1] * branch_r
    probability = float(np.vdot(selected, selected).real)
    if probability <= 1e-14:
        raise ZeroDivisionError("postselection probability is zero or unresolved")
    reduced = np.einsum("rpa,rqa->pq", selected, selected.conj()) / probability
    reduced = (reduced + reduced.conj().T) / 2.0
    eig = np.linalg.eigvalsh(reduced)
    if eig.min(initial=0.0) < -1e-9 or not np.isclose(np.trace(reduced), 1.0, atol=1e-9):
        raise ArithmeticError("reduced probe state is not a density matrix")

    # Conditional recoil/probe means from the unnormalised selected pure state.
    selected_norm = selected / np.sqrt(probability)
    recoil_after = _momentum_expectation(selected_norm, 0)
    probe_after = _momentum_expectation(selected_norm, 1)
    recoil_before = _momentum_expectation(recoil, 0)
    probe_before = _momentum_expectation(probe_state, 0)
    pointer_overlap = complex(np.vdot(pointer_initial, pointer_r))
    return StinespringResult(
        probability=probability,
        probe_state=reduced,
        joint_recoil_probe_pointer=selected,
        total_momentum_initial=initial_total,
        total_momentum_conditional=recoil_after + probe_after,
        source_recoil_shift=recoil_after - recoil_before,
        probe_shift=probe_after - probe_before,
        pointer_branch_overlap=pointer_overlap,
    )


def complete_stinespring_probabilities(
    source_amplitudes: np.ndarray,
    measurement_basis: np.ndarray,
    source_recoil: np.ndarray,
    probe: np.ndarray,
    impulse_units: int = 1,
    pointer_angle: float = 0.0,
) -> tuple[float, float]:
    """Return both source-postselection outcome probabilities."""
    basis = np.asarray(measurement_basis, dtype=np.complex128)
    if basis.shape != (2, 2) or not np.allclose(
        basis.conj().T @ basis, np.eye(2), atol=1e-10
    ):
        raise ValueError("measurement_basis must be a 2x2 unitary")
    probs = []
    for column in range(2):
        result = run_stinespring_recoil_model(
            source_amplitudes, basis[:, column], source_recoil, probe,
            impulse_units, pointer_angle
        )
        probs.append(result.probability)
    if not np.isclose(sum(probs), 1.0, atol=1e-9):
        raise ArithmeticError("complete source outcomes do not sum to one")
    return float(probs[0]), float(probs[1])
