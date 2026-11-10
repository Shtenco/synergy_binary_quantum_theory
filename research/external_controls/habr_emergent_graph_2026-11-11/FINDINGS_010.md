# Findings 010 — exact variational minimum and what it does (not) derive

Source: geometric_resonance/variation_principle.py, blob SHA 1db5bddc0f3381048705774bd6101c95c317bdc9.

## 1. Rewrite and differentiate the functional

The script defines N_c=N^(2/3), k=xN_c, K_eff=K+λk and:

    F(x)=εx+γx²−Tx ln(1/x)+β/ln(K+λxN_c).

Since −x ln(1/x)=x ln x, let A=λN_c and u=K+Ax. Then:

    F'(x)=ε+2γx+T(ln x+1)−βA/[u(ln u)²],

    F''(x)=2γ+T/x+βA²(ln u+2)/[u²(ln u)³].

For the script's positive parameters K=6, λ=10, β=3, γ=2, T=0.4, and x>0, u>6 so ln u+2>0. Every term in F'' is positive. Therefore the functional is strictly convex on x>0. Also F'(x)→−∞ as x→0+ and F'(x)→+∞ as x→∞, so there is exactly one stationary point, which is the unique global minimum.

This is a valid mathematical property of the chosen functional; the grid search is not required to establish uniqueness.

## 2. Independent continuous minimum

Solving F'(x)=0 by bisection for N=4.2×10^121 and the exact script coefficients gives:

- x*=0.0239770591844;
- F(x*)=0.00553931585893;
- F'(x*)≈4.5×10^(−15);
- F''(x*)=20.8362797808;
- k*=x*N^(2/3)≈2.895×10^79;
- p*=k*/N≈6.904×10^(−43);
- μ*=ln(k*)/lnN=0.6533451669.

The 6000-point logarithmic grid in the source gives x≈0.0239550 and curvature≈20.8520, within about 0.1% of the continuous root/curvature. The numerical minimization itself is therefore reproducible.

## 3. What the minimization does not derive

- μ_c=2/3 is assigned before the minimization through N_c=N^(2/3). The minimization does not discover the exponent 2/3.
- The exact finite-N exponent μ*=0.653345, not 2/3. The difference is ln(x*)/lnN.
- The separate article scaling p=1/(K N^(1/3)) would require pN^(1/3)=1/K=0.1666667. The variational script gives pN^(1/3)=x*=0.0239771, so Kx*=0.143862; the coefficients disagree by a factor about 6.95.
- ε, γ, T, β and λ are hard-coded in the script and are not derived there from the graph probability law or microscopic action.
- k*=xN^(2/3) is enormous. The script adds λk directly to K_eff and then uses lnN/lnK_eff as a path-length estimate. If k is total shortcut count, it cannot be added directly to a per-vertex degree; if it is mean shortcut degree, the N^(2/3) scaling conflicts with the sparse constant-shortcut-degree ensemble tested in Findings 009.
- The plot labels the minimum “critical saddle”; in this one-dimensional strictly convex functional it is a minimum, not a saddle.

## 4. Verdict

- Unique minimum of the specified F(x): PROVED analytically.
- Numerical grid minimum: REPRODUCED to about 0.1%.
- Exponent 2/3 derived by minimization: NO; inserted in N_c.
- Agreement with the article's p coefficient: FAIL.
- Physical derivation of free-energy coefficients and k interpretation: NOT ESTABLISHED.

The variational script contains a mathematically sound convex minimization, but its physical interpretation and its compatibility with the graph ensemble remain open.