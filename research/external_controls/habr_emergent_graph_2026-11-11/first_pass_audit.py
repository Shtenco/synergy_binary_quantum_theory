#!/usr/bin/env python3
"""Reproducible first-pass numerical audit. Standard library only."""
import math

K = 6.0
D = 3.0
PI = math.pi
LN_N = (K - math.log(K)) / (1.0 / D - 1.0 / PI)
LN_N_ZETA = K ** (1.5 + PI**2 / 6.0)
LN_N_PHYS = LN_N - PI * K / LN_N
DELTA_ACTUAL = LN_N_ZETA - LN_N_PHYS
DELTA_CLAIMED = PI**6 / (math.sqrt(6.0) * LN_N**2)


def lattice_partial(radius):
    total = 0.0
    for i in range(-radius, radius + 1):
        for j in range(-radius, radius + 1):
            for k in range(-radius, radius + 1):
                norm2 = i*i + j*j + k*k
                if norm2:
                    total += 1.0 / norm2
    return total


def free_energy_derivative(x, n, k=6.0, eps=1.0, gamma=2.0,
                           temp=0.4, beta=3.0, lam=10.0):
    nc = n ** (2.0 / 3.0)
    keff = k + lam * x * nc
    return (eps + 2.0 * gamma * x + temp * (math.log(x) + 1.0)
            - beta * lam * nc / (keff * math.log(keff)**2))


def bisect_root(fn, lo, hi, iterations=300):
    flo, fhi = fn(lo), fn(hi)
    if flo * fhi >= 0:
        raise ValueError("root is not bracketed")
    for _ in range(iterations):
        mid = (lo + hi) / 2.0
        fm = fn(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return (lo + hi) / 2.0


def main():
    print("RESonance")
    print(f"lnN*={LN_N:.14f}")
    print(f"N*={math.exp(LN_N):.14e}")
    print(f"lnN_zeta={LN_N_ZETA:.14f}")
    print(f"delta(zeta - phys)={DELTA_ACTUAL:.14f}")
    print(f"delta_claimed={DELTA_CLAIMED:.14f}")
    print(f"correction_relative_error_pct={(DELTA_CLAIMED-DELTA_ACTUAL)/DELTA_ACTUAL*100:.8f}")
    print()
    print("3D inverse-square lattice partial sums")
    for r in (1, 2, 4, 8, 12):
        print(f"R={r:2d} S_R={lattice_partial(r):.12f}")
    print()
    print("nonlocal zeta")
    z2, z4 = math.pi**2/6.0, math.pi**4/90.0
    print(f"zeta2={z2:.12f} zeta4={z4:.12f} ratio={z4/z2:.12f}")
    print()
    print("variation functional")
    n = 4.2e121
    x = bisect_root(lambda xx: free_energy_derivative(xx, n), 1e-8, 100.0)
    k_opt = x * n**(2.0/3.0)
    p_opt = k_opt / n
    mu_opt = math.log(k_opt) / math.log(n)
    p_article = 1.0/(6.0*n**(1.0/3.0))
    print(f"x*={x:.14g}")
    print(f"p*={p_opt:.14e}")
    print(f"mu*={mu_opt:.14f}")
    print(f"p*N^(1/3)={p_opt*n**(1.0/3.0):.14g}")
    print(f"article_p={p_article:.14e}")
    print(f"pstar/article_p={p_opt/p_article:.12f}")


if __name__ == "__main__":
    main()
