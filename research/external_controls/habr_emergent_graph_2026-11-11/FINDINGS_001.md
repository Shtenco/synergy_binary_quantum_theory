# Findings 001 — first independent recomputation

Date: 2026-11-11  
Target article: https://habr.com/ru/articles/1046730/  
Target code revision identifiers (blob SHAs):
- `physics_constants/stress_test.py`: `b617ac8f42bec5e7916dfefce576665bb68e81fc`
- `geometric_resonance/variation_principle.py`: `1db5bddc0f3381048705774bd6101c95c317bdc9`
- `matrix_basis/matrix_analysis.py`: `e58dd577e804b07d07e6c0df1a01041639ad9100`
- `spectral_mass_theory/mass_relation.json`: `6b3ca658afed9211bcbcc6c5d4e18be25e674500`
- `draft/ultimate_stresstest.py`: `154bcd4b3514f91c034328913334531773a338cc`

## 1. Geometric resonance: arithmetic reproduces, physical derivation does not

The article defines
[
R(N)=\pi\left(\frac13-\frac{K-\ln K}{\ln N}\right)
]
and imposes (R=1). Therefore, algebraically,
[
\ln N_* = \frac{K-\ln K}{1/3-1/\pi}.
]
At (K=6):
- (ln N_*=280.1115142805333);
- (N_*=4.4759473527207435\times10^{121}).

This arithmetic is correct. However, it is not yet a prediction from the graph: the order parameter (R) and critical condition (R=1) are introduced by the article. Until a microscopic graph action/dynamics derives both, this is a chosen criticality ansatz.

The article's quadratic correction
[
x^2-G_0x+\pi K=0,quad G_0=\frac{K-\ln K}{1/3-1/\pi}
]
is not algebraically equivalent to the original resonance equation. Its large root is approximately (280.0442212384), while the exact solution of (R=1) is (280.1115142805). The small root is approximately (0.06729). The quadratic correction therefore needs an independent derivation; it is not a rearrangement of (R=1).

## 2. Spectral-zeta derivation: a divergent sum is treated as finite

The article writes for 3D lattice modes:
[
\zeta_G(1)=\left(\frac{L}{2\pi}\right)^2
\sum_{n\in\mathbb Z^3\setminus\{0\}}\frac1{|n|^2}
]
and then assigns the infinite lattice sum the value (K\zeta(2)).

This is false as written. In a 3D shell of radius (r), the number of lattice points is (O(r^2dr)), while each summand is (O(r^{-2})), so the shell contributes (O(dr)). The sum diverges linearly with the cutoff.

Independent finite-cube partial sums:
- radius 1: 14.6666666667;
- radius 2: 29.8000000000;
- radius 4: 60.3459976298;
- radius 8: 121.6487439852;
- radius 12: 183.0090629761.

The steady growth is consistent with divergence, not convergence to (6\zeta(2)).

A finite graph's usual nonzero-spectrum zeta is
[
\zeta_G(s)=\sum_{i=1}^{N-1}\lambda_i^{-s}.
]
For a connected graph, at (s=0) this gives exactly (N-1), since each of the (N-1) positive eigenvalues contributes 1. It is not the article's asserted (N(2K)^{3/2}) expression. A separate regularized density-of-states quantity could be defined, but it must not be called the same zeta value without a precise definition and derivation.

## 3. The nonlocal zeta contribution changes exponent mid-derivation

The article defines (p_n=p/n^2), then writes a summand proportional to
[
p_n\ln K\frac1{n^2}=p\ln K\frac1{n^4}.
]
The sum is therefore (p\ln K\zeta(4)), not (p\ln K\zeta(2)). The subsequent replacement by (zeta(2)) is not a valid “leading order in p” step: the series index and exponent do not depend on p. Numerically:
- (zeta(2)=1.64493406685);
- (zeta(4)=1.08232323371);
- (zeta(4)/zeta(2)=0.65797362674).

This is a direct internal algebraic inconsistency.

## 4. Variation script: the exponent is inserted, and its coefficient conflicts with the article's p

The script sets (mu_c=2/3), then defines (N_c=N^{2/3}) and (k=xN_c). Consequently:
[
p=k/N=xN^{-1/3}.
]
Thus the (N^{-1/3}) exponent is built into the parameterization; the minimizer cannot independently discover that exponent.

For the exact script inputs (N=4.2\times10^{121}), (K=6), (epsilon=1), (gamma=2), (T=0.4), (eta=3), (lambda=10), the independent continuous stationarity equation is
[
F'(x)=\epsilon+2\gamma x+T(\ln x+1)
-\frac{\beta\lambda N^{2/3}}
{(K+\lambda xN^{2/3})[\ln(K+\lambda xN^{2/3})]^2}=0.
]
Bisection gives:
- (x_*=0.0239770591844);
- (p_*=6.89783526823\times10^{-43});
- (mu_* = \ln k_*/\ln N=0.653345166919).

The article's separate formula (p=1/(K N^{1/3})) implies
[
pN^{1/3}=1/K=0.1666666667,
]
whereas the variation script gives (p_*N^{1/3}=x_*=0.02397705918). Their ratio is (Kx_*=0.1438623551), i.e. a factor of about 6.95 difference. The variation script therefore does not numerically reproduce the stated coefficient (1/K).

The grid scan also depends on a finite interval and 6000 log-spaced samples; it should be compared with the continuous derivative root and convergence under grid/bounds refinement. The functional coefficients are hard-coded, not derived from the graph ensemble in this script.

## 5. The claimed NNLO difference is numerically close, but the derivation is absent

Using the article's own geometric root:
[
\ln N_{math}=280.1115142805333,
quad
\ln N_{phys}=\ln N_{math}-\frac{\pi K}{\ln N_{math}}
=280.04422123843915.
]
The claimed zeta expression evaluates to
[
\ln N_\zeta=K^{3/2+\zeta(2)}=280.0492258637182.
]
Thus the actual difference is (0.0050046252791).

The article's proposed correction (pi^6/(\sqrt6(\ln N)^2)) gives (0.0050022073305), a relative discrepancy of about (0.0483\%), not (0.0003\%) for these exact inputs. The numerical proximity is real, but the formula for the correction is not derived by the preceding spectral argument, which already contains the divergent sum and the (zeta(4)\to\zeta(2)) substitution.

## 6. Source-code / article-table mismatch: strange-quark mass

The article's displayed table gives (m_s=9.64\times10^{-30}\,kg) and formula
[
m_s=\frac{K(\ln N)^5}{\sqrt\pi N^{1/3}}.
]
The code's `constants` dictionary instead uses (1.66\times10^{-28}\,kg), and its formula implementation uses the same expression. At the stated (N\approx4.2\times10^{121}), the formula is of order (1.66\times10^{-28}\,kg), agreeing with the code dictionary, not the article table. This is a factor (sim17.2) discrepancy between article presentation and the code's reference value. It must be resolved before any aggregate accuracy statistic is trusted.

## 7. Dimensional issue: numerical formulae are not yet physical-unit derivations

The displayed right-hand sides are composed of dimensionless (K,N,\pi,\ln K,\ln N). They are dimensionless. Yet they are labelled as metres, seconds, kilograms, joules, coulombs, etc. A unit conversion scale with the appropriate dimensions must be derived and specified. Without it, the equations are numerical recipes tied to an implicit unit convention, not dimensionally complete physical equations. This issue applies to the dimensional constants and particle masses/lifetimes; it does not apply to dimensionless ratios.

## 8. Current verdict

- Arithmetic of the geometric resonance root: **PASS**.
- Physical derivation of the resonance condition: **NOT ESTABLISHED**.
- 3D inverse-square lattice sum as finite: **REFUTED AS WRITTEN**.
- Nonlocal series (zeta(4)\to\zeta(2)): **REFUTED AS WRITTEN**.
- Variation script independently deriving exponent (2/3): **NO — exponent is hard-coded**.
- Variation script reproducing (p=1/(K N^{1/3})): **FAIL**.
- Claimed 80+ constants accuracy: **NOT YET REPRODUCED; source/table mismatch already found**.
- Fundamental physical interpretation: **OPEN**.

No conclusion about the entire research direction follows from these findings alone. They invalidate specific derivations/claims as currently written and require correction or replacement.
