# Findings 003 — complete mass-ratio JSON audit

Source artifact: spectral_mass_theory/mass_relation.json
Blob SHA: 6b3ca658afed9211bcbcc6c5d4e18be25e674500

## 1. Particle count and relation count

The JSON contains 20 unique particle records and 380 directed pairwise relations. Since 20×19=380, every ordered pair of distinct particles is included. All 190 unordered pairs have a reciprocal entry; none are missing.

Therefore 380 rows are not 380 independent confirmations. They are the complete directed ratio matrix generated from 20 particle masses. At most 19 independent dimensionless mass ratios exist for 20 nonzero masses after choosing one reference mass. Reciprocal entries are algebraically dependent.

The article says that exact algebraic relations were established between 32 particles (article line 1183), but this JSON snapshot contains 20 particle records. This is a source/data-version mismatch that must be resolved.

## 2. Independent reconstruction from particle mass records

For all 380 rows, the JSON experimental_value is exactly reconstructed by dividing the two corresponding particle.mass entries. The maximum floating-point reconstruction error is at numerical roundoff. This confirms that the 380 experimental ratios are derived from the 20 mass records, not 380 additional experimental inputs.

## 3. The teta field is inverted relative to the article's own definition

The article states that teta = experimental_value / model_value (article line 1213). In the JSON, all 380 stored teta values instead equal model_value / experimental_value. The first row makes this explicit:

- model_value = 0.28209479177387814;
- experimental_value = 0.28267539990381657;
- stored teta = 0.9979460252638327 = model/experimental;
- experimental/model = approximately 1.00205820223.

All 380 rows match model/experimental and all 380 disagree with the stated experimental/model convention. This reverses the sign/direction of the reported deviation, although closeness to 1 is numerically preserved under reciprocation.

## 4. Formula-string parser audit

All 380 formula strings were independently parsed using ordinary arithmetic precedence with K=6 and pi=π. There are 338 distinct formula strings.

- 378/380 strings parse successfully.
- 2 strings are malformed and fail parsing because of unmatched/incorrect parentheses.
- 37/380 parsed rows do not reproduce the stored model_value under the literal formula string (relative mismatch greater than 1e-10).

Example: the JSON string for π-meson / Ω⁻-baryon is:

    K**(3/2) / 4 * pi**3 * 2**(1/2)

Under ordinary operator precedence this evaluates to approximately 161.11334, while the stored model_value is 0.08379194439. The stored value corresponds instead to:

    K**(3/2) / (4 * pi**3 * 2**(1/2))

Thus the formula serialization omits denominator parentheses. Similar patterns occur elsewhere. A human may infer the intended grouping from context, but the machine-readable formula is not reproducible as written. Formula strings must be corrected and versioned; the audit must not silently repair them.

Two malformed strings include examples with unbalanced parentheses in the proton/Z-boson and DT/Omega-baryon rows. They must be fixed at source, not patched only in the audit.

## 5. Statistical interpretation

The 380 rows consist of 190 reciprocal pairs. Their errors are also not independent because all ratios share the same 20 particle mass inputs and many share the same fitted basis K, pi, sqrt(2), sqrt(3). Counting all 380 as separate successes exaggerates the effective sample size.

The proper next test is particle-level: freeze the formula-selection rule and predict masses for held-out particle species or independent observables, without fitting the mass basis or N on those held-out values. A pairwise-ratio table generated after choosing 20 masses cannot substitute for that test.

## 6. Current verdict

- JSON particle records: 20, not the article's stated 32.
- Pairwise relation rows: 380, but only 190 reciprocal unordered pairs and at most 19 independent ratios.
- Experimental ratios: exactly derived from the 20 particle mass records.
- teta definition: inverted in JSON relative to article prose, 380/380 rows.
- Formula serialization: 37 literal-expression mismatches and 2 malformed expressions.
- Pairwise mass relation evidence as 380 independent confirmations: REJECTED.