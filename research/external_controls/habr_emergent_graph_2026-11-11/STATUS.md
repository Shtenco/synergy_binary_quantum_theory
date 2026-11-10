# Audit status — 2026-11-11

## Branch isolation

- Repository: Shtenco/synergy_binary_quantum_theory
- Branch: audit/habr-emergent-graph-independent-control-20261111
- Base: main
- Canonical BQG branch/main was not modified.

## Artifacts committed

- README.md — independence contract and audit sequence.
- FORMULA_LEDGER.md — claim/formula checklist.
- CONSTANT_FORMULA_REGISTER.md — 77 individually numbered displayed constant/mass/lifetime/ratio formulas transcribed from the article.
- independent_recompute.py and first_pass_audit.py — dependency-free numerical checks.
- audit_mass_relations.py — safe AST evaluator and consistency checker for the mass-ratio JSON.
- FINDINGS_001.md — geometric resonance, divergent 3D sum, zeta(4)/zeta(2), variation scaling, units.
- FINDINGS_002.md — N provenance conflicts, constants/code mismatch, strange-quark mismatch, fitting protocol.
- FINDINGS_003.md — all 380 JSON mass-ratio rows, teta inversion, formula-string parser audit, dependence structure.
- FINDINGS_004.md — neutron lifetime comparison against PDG 2026.
- FINDINGS_005.md — finite graph zeta identities, TSCO, current, unitarity, metric, mass scaling.

## Confirmed initial quantitative checks

- Geometric resonance at K=6: lnN=280.1115142805333; N=4.4759473527207435e121.
- Article's 3D inverse-square lattice sum diverges; finite cube partial sums grow from 14.6667 at radius 1 to 183.0091 at radius 12.
- Nonlocal series shown as Σ n^-4 is ζ(4), not ζ(2).
- Variation script has μc=2/3 hard-coded and returns x*=0.0239770591844; its p coefficient differs from the separate article value 1/K by factor Kx*=0.1438623551.
- stress_test.py has 82 reference constants and 83 formula outputs; B_s_meson is unmatched.
- The article, stress-test script, variation script and JSON use different N values.
- The displayed strange-quark formula evaluates to about 1.678e-28 at JSON N, while the article table prints 9.64e-30.
- mass_relation.json contains 20 particles and 380 directed relations (190 reciprocal pairs), not 32 particles as the article says.
- All 380 teta values are model/experimental, opposite to the article prose definition experimental/model.
- Literal standard-precedence parsing of JSON formula strings: 378 parse, 37 do not reproduce stored model_value, 2 fail parsing.
- Current PDG 2026 neutron mean life is 878.3±0.4 s; the article's prediction is around 879.8–880.0 s depending on N convention.

## Work still open — do not declare final verdict

1. Reproduce all 82 constants exactly by executing a frozen copy of the source script at each declared N; compare article, JSON, and code outputs row by row.
2. Independently source and validate every experimental reference value, uncertainty, renormalization scheme, and measurement method.
3. Rebuild the 55×6 exponent matrix from the original source; recompute exact rank, Smith normal form, integer kernel, LLL output, and all 29/49 claimed identities.
4. Recompute every one of the 77 displayed formulas at every N convention, and check dimensions before comparing numerical values.
5. Reproduce the variation functional's grid scan, continuous stationary point, second derivative, parameter sensitivity, and alternative exponents.
6. Specify and simulate the exact random graph ensemble; compute finite-graph Laplacian spectra, heat traces, zeta functions, spectral dimensions, and graph distances over size sweeps.
7. Prove or refute TSCO existence/uniqueness, local continuity, unitarity, continuum limit, metric signature, and particle interpretation.
8. Only then compare the surviving mechanism with BQG's binary-edge / Regge / causal-slab construction.

All findings are provisional, claim-specific, and reproducible from the committed scripts and source hashes. No complete-theory verdict has been issued.