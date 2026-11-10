#!/usr/bin/env python3
"""Independent analytic check of variation_principle.py."""
import math

N = 4.2e121
K = 6.0
EPS = 1.0
GAMMA = 2.0
TEMP = 0.4
BETA = 3.0
LAMBDA = 10.0
NC = N ** (2.0 / 3.0)
A = LAMBDA * NC


def free_energy(x):
    u = K + A*x
    return EPS*x + GAMMA*x*x + TEMP*x*math.log(x) + BETA/math.log(u)


def d1(x):
    u = K + A*x
    ell = math.log(u)
    return EPS + 2*GAMMA*x + TEMP*(math.log(x)+1) - BETA*A/(u*ell*ell)


def d2(x):
    u = K + A*x
    ell = math.log(u)
    return 2*GAMMA + TEMP/x + BETA*A*A*(ell+2)/(u*u*ell**3)


def bisect(f, lo, hi, steps=300):
    flo = f(lo)
    for _ in range(steps):
        mid = (lo+hi)/2
        fm = f(mid)
        if flo*fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return (lo+hi)/2


def main():
    x = bisect(d1, 1e-8, 100.0)
    k = x*NC
    p = k/N
    mu = math.log(k)/math.log(N)
    p_article = 1/(K*N**(1/3))
    print(f"x_star={x:.14g}")
    print(f"F_star={free_energy(x):.14g}")
    print(f"F_prime={d1(x):.4e}")
    print(f"F_second={d2(x):.12g}")
    print(f"k_star={k:.12e}")
    print(f"p_star={p:.12e}")
    print(f"mu_star={mu:.12g}")
    print(f"mu_c_input={2/3:.12g}")
    print(f"p_star_over_article_p={p/p_article:.12g}")
    print("For x>0, F''(x)>0 because all coefficients are positive and")
    print("ln(K+A*x)+2>0. Therefore this specified functional is strictly convex")
    print("and has at most one stationary point; dF/dx changes from -inf to +inf,")
    print("so the stationary point is the unique global minimum.")


if __name__ == "__main__":
    main()
