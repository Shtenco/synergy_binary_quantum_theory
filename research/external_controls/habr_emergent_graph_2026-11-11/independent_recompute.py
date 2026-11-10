#!/usr/bin/env python3
"""Independent, dependency-free checks for the first equations in Habr article 1046730.

This is not a reproduction of the author's code. It deliberately implements only
closed-form arithmetic and finite partial sums, with no fitting to experimental data.
Run: python independent_recompute.py
"""
from __future__ import annotations
import math

K = 6
D = 3
ZETA2 = math.pi**2 / 6
ZETA4 = math.pi**4 / 90


def geometric_resonance(k: float, d: float = 3.0) -> dict[str, float]:
    """Recompute R=1 solution from the article's own definition."""
    denominator = d**-1 - math.pi**-1
    if denominator == 0:
        raise ValueError("critical denominator is zero")
    ln_n = (k - math.log(k)) / denominator
    return {"lnN": ln_n, "log10N": ln_n / math.log(10),
            "log10N_over_3": ln_n / (3 * math.log(10))}


def resonance_residual(k: float, ln_n: float, d: float = 3.0) -> float:
    return math.pi * (1.0 / d - (k - math.log(k)) / ln_n) - 1.0


def cubic_lattice_inverse_square_partial(R: int) -> float:
    """Sum 1/|n|^2 over nonzero integer points in [-R,R]^3."""
    total = 0.0
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            for k in range(-R, R + 1):
                q = i*i + j*j + k*k
                if q:
                    total += 1.0 / q
    return total


def main() -> None:
    r = geometric_resonance(K, D)
    print("GEOMETRIC RESONANCE")
    print(f"K={K}, d={D}")
    print(f"ln N* = {r['lnN']:.12f}")
    print(f"log10 N* = {r['log10N']:.12f}")
    print(f"N* = exp(lnN*) = {math.exp(r['lnN']):.12e}")
    print(f"R-1 residual = {resonance_residual(K, r['lnN'], D):+.3e}")
    print()

    claimed_ln_n_zeta = K ** (1.5 + ZETA2)
    print("ARTICLE'S ZETA EXPRESSION (ARITHMETIC ONLY)")
    print(f"zeta(2) = {ZETA2:.12f}")
    print(f"K^(3/2+zeta(2)) = {claimed_ln_n_zeta:.12f}")
    print(f"delta vs resonance lnN = {claimed_ln_n_zeta-r['lnN']:+.12f}")
    print("This verifies only numerical evaluation, not the spectral derivation.")
    print()

    print("3D LATTICE SUM: FINITE CUBE PARTIAL SUMS")
    previous = None
    for radius in (1, 2, 4, 8, 12):
        value = cubic_lattice_inverse_square_partial(radius)
        increment = float("nan") if previous is None else value - previous
        print(f"R={radius:2d}: S_R={value:.12f}, increment={increment:.12f}")
        previous = value
    print("Shell estimate: number of lattice points at radius r is O(r^2);")
    print("each term is O(r^-2), hence shell contribution is O(1) per dr.")
    print("Therefore the 3D sum over all nonzero lattice points diverges.")
    print()

    print("FINITE GRAPH ZETA AT s=0")
    print("For a connected finite graph with N vertices, there are exactly N-1")
    print("positive Laplacian eigenvalues. Thus zeta_G(0)=sum_{i=1}^{N-1}1=N-1.")
    print("This is not asymptotic to N*(2K)^(3/2) with fixed K as an identity.")
    print()

    print("NONLOCAL SERIES")
    print(f"zeta(2) = {ZETA2:.12f}")
    print(f"zeta(4) = {ZETA4:.12f}")
    print(f"zeta(4)/zeta(2) = {ZETA4/ZETA2:.12f}")
    print("The displayed series sum n^-4 is zeta(4), not zeta(2).")


if __name__ == "__main__":
    main()
