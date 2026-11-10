# Findings 014 — calibration, NLO scoring, and chi-square are not independent validation

Source: physics_constants/constant_matrix_analysis.py, blob SHA 79b5e570ab990c6143ace748881ba20203de049c.

## 1. The reported R² and AIC calculations square an already-squared residual

In extract_nlo_corrections(), NumPy's:

    coef, res, _, _ = np.linalg.lstsq(X, y, rcond=None)

returns `res` as the residual sum of squares (SSE), not the vector of residuals. The code then calculates:

    R² = 1 − sum(res**2) / TSS
    AIC = n*log(sum(res**2)/n) + penalty

That squares SSE a second time. Correct expressions (when `res` is populated) are:

    R² = 1 − SSE/TSS
    AIC = n*log(SSE/n) + 2k + constant

Alternatively, calculate the residual vector `y-X@coef` and square its entries exactly once. The present code's R² and AIC values are therefore mathematically incorrect and can select the wrong NLO model.

## 2. NLO performance is measured in-sample after a prior full-sample calibration

logarithmic_calibration_and_deviations() first fits slope and intercept by regressing all ln(CODATA) values on all ln(ETI) values. It then constructs residual deviations for the same observables. extract_nlo_corrections() fits its NLO/NNLO regressions on those same residuals and reports R²/AIC on the same rows.

This is two stages of fitting on the same targets, with no held-out set. It measures descriptive in-sample fit, not out-of-sample prediction. The intercept and slope also absorb overall dimensional scale and exponent mismatch, so a high R² cannot validate the uncalibrated physical formulas.

## 3. The regression mixes dimensionful SI values with a dimensionless monomial basis

compute_predicted_value() evaluates products of dimensionless primitives (sqrt2, sqrt3, pi, lnK, lnN, N^(1/3)). These outputs are dimensionless numbers. The code then takes logarithms of CODATA values expressed in different SI units (kg, s, J, C, m, etc.) and regresses them against the dimensionless outputs.

A logarithm of a dimensionful quantity is not invariant: changing kg to g or seconds to milliseconds adds a constant to the log value. The fitted intercept therefore changes with unit conventions. The regression is not a dimensionally defined physical model unless every quantity is first divided by a specified reference unit/scale and that scale is itself part of the theory.

## 4. The chi-square uncertainty model assumes independence

compute_identity_value() estimates the log variance as:

    Var(log I) = Σ_i (v_i σ_i)^2

which assumes all input measurements are statistically independent. But many entries are correlated or definitionally linked: h=2πℏ; Planck units depend on G, c and ℏ; ε₀, μ₀, Z₀ and α are related by exact electromagnetic identities; and particle ratios are derived from the same masses.

The code does not use the covariance matrix of the CODATA/PDG adjustment or propagate shared inputs. It can therefore misstate the uncertainty of an identity and its Z-score. Treating exact-defined constants as zero-uncertainty independent observations further makes some identity tests degenerate.

## 5. The classification is not a statistical test of theory predictions

chi_squared_test() classifies relations according to whether their exponent vector contains selected mass/coupling names and a coefficient-norm threshold. That classification is a heuristic label, not a pre-registered distinction between derived identities and independent physical predictions. It does not account for the fact that all kernel identities are generated algebraically from the same chosen exponent matrix.

## 6. Required replacement protocol

1. Freeze the exponent matrix and derive a saturated integer-kernel basis.
2. Freeze N and all global scales using a derivation that does not use the target constants being tested.
3. Use dimensionless observables or explicitly defined reference units with propagated uncertainty.
4. Use the full covariance matrix, including exact algebraic dependencies and shared measurement inputs.
5. Pre-register a train/test split by independent physical observable families, not random rows.
6. Fit any NLO coefficients only on training families; score on untouched held-out families.
7. Compute SSE directly from residual vectors and use correctly defined R²/AIC/BIC.

## 7. Verdict

- NLO R² formula: FAIL (SSE squared a second time).
- NLO AIC formula: FAIL (SSE squared a second time).
- NLO validation: in-sample after full-sample calibration; not predictive evidence.
- Regression across raw SI units: dimensionally ill-defined without reference scales.
- Chi-square uncertainty propagation: incomplete; covariance ignored.
- Statistical claim of predictive accuracy: NOT ESTABLISHED by this toolkit.