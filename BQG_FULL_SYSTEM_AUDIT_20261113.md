# BQG full-system audit checkpoint — 2026-11-13

## Executive result

The repository contains a mathematically proved **reduced-model propagation cone** and a growing causal-response reference/testing layer. The reduced cone is real work and is not being dismissed. It is not yet the physical Lorentz cone of the full graph-changing BQG theory.

Current machine-checked status remains:

- `structural_candidate_closed = false`
- `physicalization_complete = false`
- `experimentally_confirmed = false`
- all 9 required physicalization gates and all 4 supplemental physical gates remain `open_physical`.

No physical gate was closed by this audit.

## 1. Proven causal-cone result

The reduced homogeneous free-fermion / staggered-volume model has dispersion

\[
E_\pm(k)=\pm\sqrt{4m^2+9\cos^2 k}.
\]

Its sharp asymptotic ballistic front is

\[
v_{\mathrm{front}}(m)=\sqrt{4m^2+9}-2|m|.
\]

The nearest complex-momentum branch point gives

\[
\mu_c=\operatorname{arsinh}(2|m|/3)=\xi^{-1},
\]

and the static/dynamic identity is

\[
v_{\mathrm{front}}/3=e^{-1/\xi}.
\]

At criticality, \(m=0\), the exact propagator is \(U_r(t)=i^rJ_r(3t)\), with front speed 3 in the model's own units. For every fixed positive contour displacement inside the analytic strip, the propagator has an exponential exterior bound; taking the contour displacement to zero recovers the exact group-velocity front.

Local execution of `scripts/bqg_sharp_analytic_causal_cone_gate.py` passed 25 contour-grid comparisons: five masses (0.2, 0.5, 1, 2, 5) and five analytic-strip fractions per mass, plus the massless critical limit and the static/dynamic identity.

**Scope:** this proves the cone in the stated reduced free-fermion model. It does not prove the continuum Lorentz cone of full physical BQG, exact microcausality for the interacting graph-changing theory, or a measured speed of light.

## 2. Independent LeMaitre–Perche comparison: limited but useful

The pinned external 1+1 reconstruction metric has null slopes

- \(-0.9999391238985609\)
- \(+0.9999389681654915\).

The reduced BQG massless front is 3 in its own units; dividing by 3 gives a dimensionless value of 1. The reference comparison passes the declared tolerance.

This is a **dimensionless consistency check only**. The external work uses a massless scalar field on an already prescribed Minkowski background, supplies physical coordinates to the agent, and reconstructs a 1+1 metric from a QMI-based objective. The stored metric was obtained using a neighbour-ridge mask and an unweighted quadratic fit with \(g_{xx}=1\). The offline gate does not rerun the external notebook or five-agent training. It therefore does not establish common microscopic dynamics, BQG geometrogenesis, photon propagation, or physical unit calibration.

## 3. Physical history / causal response layer

The W5 exact Hermitian spectral certificate also passes. It proves the finite Euclidean Gram lower bound \(\\lambda_{\\min}(G_{EE})>32.41356\), hence \(\\sigma_{\\min}(A)>5.6932907883\). Therefore a cancellation in \(H=-(2/3)A+Q\) requires \(\\|Q\\|_{op}>3.7955271922\), while a cancellation in \(H=A+\\lambda B\) requires \(|\\lambda|\\|B\\|_{op}>5.6932907883\). The Lorentzian map \(B\) is not present in this certificate; these are necessary thresholds, not a measured cancellation or physical joint-kernel result.

## 3. Physical history / causal response layer

Additional local controls were executed after the first draft of this checkpoint:

- The finite master-projector control passed with Hilbert dimension 8, four constraints, a two-dimensional common kernel, common-kernel residual \(1.35\times10^{-15}\), and projector agreement across four positive metric choices at approximately \(10^{-15}\).
- The C8 relational-history positive control passed for the externally supplied clock shift and \(R=J\); all eight conditional states matched \(R^t|\psi_0\\rangle\).
- The relational metric-source positive control passed: \(Z(j_x,j_z)=\\cosh\\sqrt{j_x^2+j_z^2}\), with a rank-2 finite shape Hessian.
- The projector-only causality NO-GO passed: two positive master operators have the same full physical projector but different support-graph distances (1 versus 2). Therefore the physical projector alone cannot determine causal order.
- The graph-changing causal-order theorem gate passed 1,392 master-sector-pair checks, 1,392 generic lower-bound saturations, and 120 relational-locality zero checks. Under its declared one-move locality assumption it establishes \(d_M\\ge\lceil d_C/2\rceil\) and master-flow locality.

These are valuable exact finite structural controls. The master heat-kernel parameter is still not physical spacetime time; the C8 clock and \(R=J\) are external positive controls; the projector no-go shows why the zero-projector limit alone loses the causal support information. None of these results constructs the BQG physical retarded kernel.

The physical-history reference engine now covers:

- physical-environment Gram positivity for a supplied metric;
- decoherence-functional positivity and interference controls;
- support depth versus coherent depth;
- synthetic channel-front comparisons across refinement levels.

The inputs in the reference tests are supplied/synthetic. They are not actual BQG physical states or retarded kernels. A kinematical partial trace, QMI landscape, symmetric correlator, master-constraint heat-kernel parameter, or constraint-resolvent spectral variable is not a substitute for a physical retarded response.

The 4 supplemental physical gates are present and correctly remain open:

- `PHYSICAL_ENVIRONMENT_OVERLAP`
- `PHYSICAL_BORN_MEASURE`
- `PHYSICAL_RELATIONAL_CLOCK_RESPONSE`
- `PHYSICAL_COMMON_CAUSAL_CONE`

The required 9 physical gates also remain open: physical projector/history, connected inter-block history, physical TT kernel, IR six-vector, common scale calibration, dynamical Maxwell kernel, physical background cosmology, physical scalar cosmology, and lensing dynamics closure.

## 4. Repository-wide truth-ledger audit

Local execution on the staged source snapshot:

- `scripts/verify_theory_gates.py`: valid; 55 structural gates — 23 `proved`, 21 `tested_finite`, 4 `conditional`, 4 `external_extension`, 2 `experimental_test`, 1 `no_go`. It correctly reports `structural_candidate_closed=false`.
- `scripts/verify_physicalization_gates.py`: valid; 23 physicalization ledger entries — 2 `proved`, 6 `tested_finite`, 13 `open_physical`, 2 `experimental_test`. `physicalization_complete=false`.
- `scripts/audit_core_scope.py`: PASS; 525 active text/code files and 221 active Python files; historical archive quarantined; no retired vocabulary or retired-module imports.

The single core `no_go` is `PURE_ANCILLA_NODE_REFINEMENT`: a fixed pure four-qubit Gauss-singlet ancilla cannot induce the canonical S4-equivariant map. This is a no-go for that specific shortcut, not a refutation of the replacement bilinear construction. The bilinear refinement gate locally passes with S4 covariance error \(9.996\times10^{-16}\), projection leakage \(3.43\times10^{-31}\), and rank 2. Whether this no-go should remain classified as a core gate or be recorded as a rejected construction route is a theory-governance decision; do not silently reclassify it.

A documentation inconsistency was found and corrected in this draft branch: the current top of `OPEN_PROBLEMS.md` declared `structural_candidate_closed=true`, while the executable verifier returned false. `OPEN_PROBLEMS.md` and the dated header of `THEORY_STATUS.md` now agree with the 2026-11-13 machine ledger; the 2026-09-29 true declaration is explicitly retained as historical rather than silently erased.

## 5. Production W5 source staging and reproducibility

The W5 staging artifact from workflow run `43442708226` was downloaded and inspected as data, not executed as an archive.

- outer artifact SHA-256: `998a23e7d931a17dd210f9500f0164ca4713eb8c73cf0f89d16a92862203a54f`
- nested source archive SHA-256: `fa388676eb8489756c693ad96ee956d19c2754fdd82e85812a5bc20a5dd9b6c5`
- source commit: `e7f6b983e2990b1c268d9164f7f2e438f54e9e6e`
- source tree: `60e32745e0fb630964714876775e8cde20e1adab`
- all 760 manifest-listed source files passed byte-count and SHA-256 verification.

The staging runner recorded Python 3.12.16 with `numpy==2.5.4`, `sympy==1.14.0`, and `mpmath==1.3.0`. The local audit container used Python 3.13.5, NumPy 2.3.5, SymPy 1.14.0, and mpmath 1.3.0; therefore local numerical execution is not bit-for-bit environment reproduction. The W5 production/L1/L2 workflows use Python 3.11 and install unpinned latest NumPy/SymPy. **Dependency/runtime pinning across staging and production workflows is a reproducibility task still open.**

Local execution on the extracted source snapshot:

- `python -m compileall -q scripts tests`: exit 0; one non-fatal Python 3.13 SyntaxWarning for an invalid `\o` escape in a docstring at `scripts/q2_relational_history_projector_gate.py:15`.
- 56 focused W5/adjoint/common-carrier/dispatch regression tests: PASS.
- analytic reduced cone gate: PASS.
- LeMaitre finite-reference comparison: PASS.
- physical-history reference engine: executed successfully.
- theory and physicalization truth validators: PASS.
- core-scope audit: PASS.

The full local V10 W5 production closure did not complete within the 300-second execution limit, even for a one-seed/one-node probe. This is recorded as `NOT RUN TO COMPLETION`, not as a scientific failure or a successful result.

The bounded GitHub Actions production closure study on `research/bqg-w5-production-audit-20261113` (run `43444229133`) completed with **all three cutoff jobs failing** at `Jmax2=10,12,14`. They failed before writing closure JSON, at the same fail-closed invariant in `restore_scalar_charge_M`: the sine-ordered `H_E` path was expected to preserve the charged endpoint's `J,M`, but at node 0 an input `J2=0, M2=0` produced an output label `(J2,M2)=(2,0)`. The stack passes through `complete_HE_sine` and the cached Wigner–Eckart reduction. This is a concrete production blocker: either the charged-sector scalar-selection assumption/cache reduction is invalid for this operator, or the finite operator/pruning/cutoff path violates the claimed scalar Ward identity. It is **not** a physical no-go theorem about BQG; it is a failed internal production invariant that must be resolved before accepting V10 closure. No closure artifact was produced.

To discriminate between a pruning-induced symmetry leak and a deeper operator/canonicalization issue, a one-node `shadow_no_prune` diagnostic was launched separately (run `43445958758`, `Jmax2=10`, source node 0, all five seed labels with a five-state/depth-zero closure cap). Its artifact must be inspected before choosing the repair. A more precise paired context probe (run `43446220696`) now runs the existing `peter_weyl_w5_scalar_failure_context_probe.py` in both `production` and `shadow_no_prune` modes for the single `K=0` seed and source node 0. It captures the exact input key, operator stage, per-sequence output-sector weights, and largest complex amplitudes while leaving the scalar guard enabled. These two artifacts will determine whether the residue appears before or only after thresholded accumulation.

## 6. Input-integrity hardening in the causal-response reference engine

Draft PR #82 adds fail-closed validation for:

- non-finite physical metrics, state vectors, class operators, history amplitudes, and response samples;
- non-finite/negative thresholds;
- fractional or boolean history depths that were previously coerced to integers.

Five regression tests cover malformed inputs. The workflow was expanded to execute both the reduced analytic cone gate and the LeMaitre finite-reference comparison in addition to the physical-history tests and fail-closed physicalization ledger.

GitHub Actions run `43443370864` completed successfully, including compile/run, analytic cone and independent reference comparison, regression tests, physicalization ledger validation, and fail-closed status enforcement. A numerical edge-case hardening then rationalized the massive front formula to avoid cancellation for large mass and defined the exact \(\mu=0\) limits; run `43444799790` also completed successfully. The expanded finite-clock/graph-causality workflow first failed in run `43445403259` because the job installed NumPy but omitted SymPy required by the C8 history control; the dependency was added and run `43445580212` then passed all steps, including master projector, relational history, relational metric source, projector-only causality no-go, graph-changing causal-order theorem, W5 reduced-cone certificate, analytic cone, LeMaitre comparison, regression tests, and fail-closed physical ledger. The non-raw `\\otimes` docstring in the C8 history script was also corrected to a raw docstring and compiled with `SyntaxWarning` treated as an error.

The separate L1 accumulator audit completed successfully in run `43444754151`; its 9 tests include 3 expected failures that characterize legacy `COMP.add` / `CK.add_state` sub-threshold nonlinearity. The L1 strict reference accumulator passes, but it is not yet wired into the production V09/V10/V11 path.

The L2 real-production-primitive audit completed successfully in run `43444729954`: 12 tests in 92.8 seconds, with 3 expected failures for threshold-boundary additivity and passing safe-domain homogeneity/additivity checks on independent basis states. These expected failures are evidence of known state-dependent threshold behavior, not green acceptance of full linearity.

The full five-node/five-seed V10 closure study remains in progress in run `43444229133` at `Jmax2=10,12,14`; its artifacts have not yet been inspected.

PR #82 is deliberately draft and unmerged.

## 7. Production W5 blockers that remain

The source audit identifies several independent issues that must not be collapsed into a single PASS:

1. Lower-level numerical pruning is instrumented only partially; the V10 result explicitly sets `lower_level_pruning_coverage_complete=false`. `shadow_no_prune` is therefore not yet a true no-pruning reference.
2. The incremental-pruning counterexample shows that deleting a running key after each sub-threshold addition is not equivalent to coherent accumulation followed by one final threshold. Its production frequency and magnitude must be measured from the actual path.
3. Spin-wall omissions, amplitude-threshold omissions, spectral zero-aware regularization, and state/depth caps are distinct error channels and need separate ledgers and convergence studies.
4. The 36x36 fixture is a two-node supplied fixture, not a five-node production W5 amplitude export.
5. A direct-adjoint consistency result is not by itself a physical Hermitian Hamiltonian, a positive physical inner product, or a physical TT kernel.

The production `ZERO_LOSS` claim remains **NOT CERTIFIED**.

## 8. Required next sequence

1. Resolve the failed charged-sector `H_E^sine` scalar-selection/Wigner–Eckart invariant. Compare the exact full-M implementation with the highest-weight cache path; do not delete or weaken the assertion until the Ward identity is independently established.
2. Inspect the one-node `shadow_no_prune` diagnostic to determine whether intermediate pruning causes the `J=0 -> J=1` residue; if it does, trace each threshold site and restore coherent accumulation before thresholding.
3. After the scalar-selection invariant is repaired, rerun the full `Jmax2=10/12/14` W5 production study and compare full complex supports, norms, spin-wall omissions, amplitude losses, closure caps, and direct-adjoint residuals.
4. Run cache-on/cache-off and production-vs-shadow comparisons on the same frozen source, spectral convention, input basis, and node ordering.
5. Complete a transitive threshold-site inventory; close the `lower_level_pruning_coverage_complete` gap before interpreting any no-pruning result.
4. Derive the physical clock or boundary-history ordering from the actual graph-changing constraints and construct the positive physical inner product/rigging map.
5. Build physical history class operators and gauge-invariant relational observables; only then compute intervention-based physical retarded kernels.
6. Compare photon, matter, scalar, vector, and TT fronts on one physical carrier and background; independently reconstruct the spectral metric and test regulator/volume/clock/state convergence.
7. Resolve the pure-ancilla no-go's role in the canonical gate policy without hiding the counterexample.

## Final audit verdict

**Reduced causal cone: proved within its reduced free-fermion model.**  
**Causal-response reference infrastructure: implemented and regression-tested.**  
**Full physical BQG causal cone: OPEN.**  
**Full physical BQG theory: not complete and not experimentally confirmed.**

This checkpoint is an evidence ledger, not a declaration of theory completion.
