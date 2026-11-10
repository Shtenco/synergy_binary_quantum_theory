# Findings 008 — the spectral-zeta route to ln N breaks at multiple algebraic steps

Article source: lines 101–121 and 147–195 in the spectral-zeta derivation.

## 1. Mellin transform includes the zero mode incorrectly

The article defines the nonzero-spectrum zeta function ζ_G(s)=Σ_{n=1}^{N−1} λ_n^(−s), with λ_0=0 excluded. But it then writes the Mellin transform using the full heat trace K(t)=Σ_{n=0}^{N−1} exp(−tλ_n). The n=0 term is exactly 1 for all t, so the integral ∫_0^∞ t^(s−1)K(t)dt diverges at infinity for Re(s)>0. The correct finite-graph identity uses the reduced heat trace K(t)−1 (or Tr' exp(−tΔ)) for the nonzero spectrum. This is a direct mismatch between the zeta definition and the heat-trace transform.

## 2. The geometric Weyl-count factor is dropped

The article derives N_modes ~ (2K)^(3/2) from λ_max=2K, then the final formula uses K^(3/2) as the geometric factor. The factor 2^(3/2) disappears without a stated normalization derivation.

At K=6:
- (2K)^(3/2)=41.56921938;
- K^(3/2)=14.69693846;
- omitted factor=2^(3/2)=2.82842712.

If the proportionality constant in the Weyl law is intended to absorb this factor, it must be derived and consistently propagated. The article moves from an asymptotic proportionality to a precise numerical equality without fixing that normalization.

## 3. The 3D local-mode sum diverges at s=1

The article writes Σ_{n∈Z^3\{0}} |n|^(−2)=Kζ(2). This sum diverges. In dimension d, the lattice sum Σ|n|^(−2s) converges only for Re(s)>d/2. Here d=3 and s=1, so 1 is below the convergence threshold 3/2. A finite cutoff or analytic continuation is needed; neither gives the asserted positive finite value Kζ(2) without a specific regularization prescription.

## 4. The nonlocal contribution loses both the ζ exponent and the small factor p

The article itself first derives:

    ζ_nonloc(1) = Σ_n p_n lnK / n²
    with p_n=p/n²
    = p lnK Σ_n n^(−4)
    = p lnK ζ(4).

That line is algebraically correct. The next line replaces ζ(4) with ζ(2) “in leading order in p”. But p is a multiplicative prefactor; taking a leading-order expansion in p does not change the summation exponent from 4 to 2. The series constants are different: ζ(4)=1.08232323371 and ζ(2)=1.64493406685.

Even if one granted the unsupported ζ(4)→ζ(2) replacement, exponentiating the contribution must give:

    exp[p ζ(2) lnK] = K^(p ζ(2)),

not K^(ζ(2)). The factor p is dropped entirely in the article's exponentiation. At N=4.2×10^121 and K=6, p=1/(K N^(1/3))≈4.7947465×10^(−42). Therefore:

- p lnK ζ(4) ≈ 9.2983×10^(−42);
- exp[p lnK ζ(4)] = 1 + 9.2983×10^(−42) + O(p²);
- K^(ζ(2)) ≈ 19.0549363.

The derived spectral multiplier is therefore approximately 1, not 19.05, under the displayed equations.

## 5. Final formula is not obtained from the preceding steps

The article concludes lnN_ζ=K^(3/2)×K^(ζ(2))=K^(3/2+ζ(2))≈280.0492. But:

- the local spectral sum used to obtain ζ(2) diverges;
- the Mellin transform includes an unremoved zero mode;
- the geometric factor changes from (2K)^(3/2) to K^(3/2) without a fixed normalization;
- the nonlocal sum gives ζ(4), not ζ(2);
- the small probability p is removed when exponentiating;
- no equation in the preceding derivation equates the resulting product to lnN.

The numerical value 280.0492 is reproducible from the final expression, but the derivation does not establish that expression. This is not a small correction: the missing p changes the nonlocal multiplier from 1+O(10^(−42)) to approximately 19.

## 6. Current verdict

- Mellin transform as printed: REFUTED for the full heat trace with zero mode included.
- 3D local sum at s=1: DIVERGENT.
- (2K)^(3/2)→K^(3/2): NORMALIZATION GAP.
- ζ(4)→ζ(2): ALGEBRAICALLY INVALID.
- exp[pζ(2)lnK]→K^ζ(2): ALGEBRAICALLY INVALID; p was dropped.
- Final lnN_ζ numerical arithmetic: PASS; spectral derivation: REFUTED AS WRITTEN.