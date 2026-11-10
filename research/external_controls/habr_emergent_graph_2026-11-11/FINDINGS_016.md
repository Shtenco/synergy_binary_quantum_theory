# Findings 016 — pairwise mass-ratio table has rank 19, not 380 independent tests

Sources:
- `spectral_mass_theory/mass_relation.json`, blob SHA `6b3ca658afed9211bcbcc6c5d4e18be25e674500`.
- `spectral_mass_theory/mass_analyzis.py`, blob SHA `957c0d8ffb576ccd09f40cb77015f449b8018e98`.
- Prior independent reconstruction: FINDINGS_003.md.

## 1. Exact dependency count

The JSON contains (n=20) particle masses and all (n(n-1)=380) ordered ratios. In logarithmic coordinates,
[
r_{ij}=\log(m_i/m_j)=\ell_i-\ell_j.
]
The complete ratio table is represented by the incidence matrix of the complete graph (K_{20}). Its rank is exactly (n-1=19): the common shift (ell_i\mapsto\ell_i+c) cancels from every ratio. Thus:
- 380 directed rows;
- 190 unordered pairs;
- 19 independent log-ratio degrees of freedom at most;
- 361 linear dependencies among the 380 directed rows.

The reciprocal relation (r_{ji}=-r_{ij}) already accounts for half the duplication; the remaining dependencies arise because ratios sharing masses obey cycle identities, e.g. (r_{ij}+r_{jk}=r_{ik}). This is exact algebra, not a statistical approximation.

## 2. Experimental ratios are reconstructed from the same mass table

FINDINGS_003 established that every `experimental_value` is the quotient of the two `particles[].mass` values. Consequently the 380 experimental ratios do not constitute an independent experimental dataset. They are a complete re-expression of the same 20 inputs.

A dependency-aware audit must therefore work at particle/observable level, retain the original measurement covariance, and only then construct ratios. A chi-square over 380 rows with diagonal errors would count the same mass measurements repeatedly and produce invalid degrees of freedom.

## 3. Missing uncertainty and provenance block a valid statistical test

The JSON stores point masses but no per-particle uncertainty, covariance matrix, source citation, measurement method, or (for quark masses) renormalization scheme and scale. It is therefore impossible to reconstruct a statistically valid covariance matrix for the ratio vector from this JSON alone.

For log-mass vector ℓ with covariance Σ, the ratio covariance must be
[
\operatorname{Cov}(B\ell)=B\Sigma B^T,
]
where each row of (B) has one (+1) and one (-1). This covariance is necessarily singular (rank at most 19), even if Σ is positive definite. Treating all 380 ratios as independent with diagonal errors is mathematically inconsistent.

## 4. Additional source-level duplication in mass_analyzis.py

The separate `mass_analyzis.py` table repeats the same experimental value for distinct labels:
- `m_DT` and `m_D0` both use (3.351077\times10^{-27}\,\mathrm{kg});
- `Ksi_0` and `Ksi_minus` both use (2.369569\times10^{-27}\,\mathrm{kg}).

These may be intentional approximations or transcription/alias errors, but they cannot be counted as independent observations without source-level justification. The file also includes derived constants and ratios alongside particle masses in one pooled error analysis; these share inputs by construction. A row-wise Shapiro, t-test, ANOVA, or mean-error interpretation assumes an independence structure that the script does not establish.

## 5. Corrected test protocol

1. Freeze one canonical list of physical observables, not every algebraic ratio.
2. Attach a source, central value, uncertainty, scheme/scale, and covariance provenance to each observable.
3. Partition calibration observables from genuinely held-out particle species or measurements before choosing formula families or (N).
4. Construct ratios only after the covariance of the primary measurements is defined.
5. Use (B\Sigma B^T), a rank-aware likelihood, or a nonredundant 19-dimensional contrast basis; do not invert the 380×380 singular covariance as if it were full rank.
6. Report the effective number of independent constraints and all formula-selection degrees of freedom.

## 6. Verdict

- “380 pairwise relations” as “380 independent confirmations”: **REFUTED**.
- Maximum independent log-ratio rank for 20 masses: **19**.
- Statistical significance from the present JSON alone: **NOT IDENTIFIABLE** (uncertainties/covariance/provenance absent).
- Duplicate experimental entries in `mass_analyzis.py`: **CONFIRMED**, interpretation unresolved.
- Particle-level blind prediction: **NOT YET PERFORMED**.

Reproducible dependency checker: `audit_mass_dependency.py`.
