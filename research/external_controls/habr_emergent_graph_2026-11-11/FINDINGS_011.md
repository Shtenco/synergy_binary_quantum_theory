# Findings 011 — larger sparse-graph control reverses the preliminary ln⁻² gap claim

## 1. Correction to Findings 009

Findings 009 used only L=8,12,16 and two realizations per size. Its table suggested that gap·(lnN)² was roughly stable. That conclusion was premature. Extending the same graph law to L=16,24,32,48 with three seeds per size shows gap·(lnN)² increasing systematically. The random-graph spectral gap itself trends toward a positive constant over this range, not toward zero as 1/(lnN)². Findings 009 must therefore be read as an exploratory small-size observation, not evidence for the claimed asymptotic law.

## 2. Larger sparse-graph experiment

To avoid the O(N²) pair scan in the first prototype, the new script uses a sparse Poissonized edge sampler. Each endpoint proposes a nonlocal edge at rate p/(2r²); an undirected edge is present if either endpoint proposes it. Therefore:

    P(edge) = 1 − exp(−p/r²) = p/r² + O((p/r²)²).

For the tested sizes the maximum relative deviation from p/r² is below about 0.3%. The six local cubic edges remain deterministic. This preserves the target kernel to controlled second-order accuracy while allowing N up to 110,592.

Implementation: sparse_scaling_control.py. Parameters: p=1/(6L), periodic L³ torus, three independent seeds at each size.

## 3. Sampled random-graph spectral gaps

| L | N | Mean degree | Mean λ₂ | SD(λ₂) | Mean λ₂·(lnN)² |
|---:|---:|---:|---:|---:|---:|
| 16 | 4,096 | 7.1317 | 0.62310 | 0.00588 | 43.11 |
| 24 | 13,824 | 7.1794 | 0.55553 | 0.01099 | 50.50 |
| 32 | 32,768 | 7.2025 | 0.53742 | 0.00217 | 58.10 |
| 48 | 110,592 | 7.2288 | 0.52353 | 0.00258 | 70.61 |

The normalized quantity rises by about 64% from L=16 to L=48. It is not approaching a stable constant over this range. Meanwhile λ₂ decreases from 0.623 to 0.524 and appears to flatten. A simple least-squares fit to the four size means gives:

    λ₂(L) ≈ 0.4649 + 2.4315/L

with residual sum of squares 2.29×10⁻⁴. By contrast, a zero-intercept c/(lnN)² fit has residual sum of squares 3.85×10⁻², about 168 times larger. This is descriptive finite-size model comparison, not a proof of the limit.

## 4. Exact spectrum of the expected Laplacian

The ensemble is translation-invariant on the torus. Its expected adjacency kernel is a convolution, so its expected combinatorial Laplacian is diagonalized exactly by discrete Fourier modes. For wavevector q:

    λ̄(q) = 2Σ_{a=1}³[1−cos(q_a)]
            + p Σ_{r∈T_L, |r|>1} [1−cos(q·r)]/|r|².

The minimum nonzero eigenvalue of this deterministic expected operator was computed by FFT over all Fourier modes, for L=8,12,16,24,32,48,64,96,128,192,256,384,512:

| L | λ₂(E[L]) | λ₂(E[L])·(lnN)² |
|---:|---:|---:|
| 16 | 0.876081 | 60.61 |
| 32 | 0.765885 | 82.79 |
| 64 | 0.737790 | 114.85 |
| 128 | 0.730718 | 154.82 |
| 256 | 0.728946 | 201.73 |
| 512 | 0.728502 | 255.16 |

The expected-operator gap converges numerically toward a positive constant. For the lowest coordinate Fourier mode q=(2π/L,0,0), the Riemann-sum limit of the nonlocal term is:

    C = (1/6) ∫_{[-1/2,1/2]^3} [1−cos(2πu_x)]/(u_x²+u_y²+u_z²) d³u.

The local-lattice contribution tends to zero as L⁻². Eight independent scrambled Sobol integrations (2²⁰ points each) give C≈0.72835357 with estimated standard error 5.2×10⁻⁷. The direct finite-L FFT values approach this coefficient (λ₂(E[L])=0.72850245 at L=512). The minimum nonzero Fourier mode was a coordinate-axis mode for all tested L. A proof that no other mode yields a smaller asymptotic value is still needed before calling this the exact λ₂ limit.

Important: λ₂(E[L]) is not E[λ₂(L)]—the latter is nonlinear—so this expected-operator calculation does not replace the sampled random-graph experiment. The independent samples point in the same qualitative direction: a positive gap plateau is more plausible than a ln⁻² collapse.

## 5. What remains true about p∝N⁻¹/³

The expected shortcut degree is:

    E[d_shortcut] = p Σ_{r≠0,local} 1/|r|².

In 3D the lattice shell count grows as r², so the sum grows linearly with L. Hence p∝L⁻¹=N⁻¹/³ keeps the expected shortcut degree O(1). The finite-torus shell sums confirm this: the extra expected degree rises slowly from about 0.96 at L=8 to about 1.25 at L=96, approaching a finite value.

This part survives conditionally. It is a consequence of the assumed three-dimensional shell geometry and r⁻² kernel; it does not independently establish either assumption.

## 6. Revised verdict

- p∝N⁻¹/³ for finite expected shortcut degree: CONDITIONAL DERIVATION SURVIVES.
- Preliminary claim λ₂∼1/(lnN)² from L≤16: RETRACTED as unsupported.
- Larger random-graph samples L=16…48: GAP PLATEAU IS MORE CONSISTENT with the data than ln⁻² decay.
- Exact expected-Laplacian spectrum through L=512: GAP→approximately 0.7284, not zero.
- Asymptotic theorem for the random graph's E[λ₂]: NOT PROVED; needs larger sizes, concentration bounds, and rigorous random-operator analysis.
- No inference from this graph experiment to Lorentzian spacetime, particle masses, or BQG is warranted.