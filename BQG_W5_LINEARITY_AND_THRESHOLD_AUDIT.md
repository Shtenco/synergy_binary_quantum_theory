# W5 linearity and amplitude-threshold audit

Date: 2026-11-28

## Executive result

The current W5 adjoint-closure pipeline is not yet an exact finite-dimensional linear-operator construction. It contains hard amplitude thresholds at several layers, not only in the v10 breadth-first closure. For a fixed positive threshold, coefficient pruning is not homogeneous:
\[
P_\tau(\lambda x)\ne \lambda P_\tau(x)
\]
for some nonzero \(x\) and \(\lambda\). It can also make sparse accumulation depend on the order/grouping of contributions. Consequently, a thresholded reachable-set closure and its matrix-adjoint residual are useful finite numerical diagnostics, but they do not certify the exact invariant subspace or the exact W5 operator.

This is a mathematical limitation of the current numerical realization, not a proof that the underlying BQG operator is nonlinear.

## Source-level findings on the pinned research branch

The following paths were inspected at branch head 625090067a8f7999b50d1ebef267b826ca16f9f9:

| Layer | Source path | Observed threshold behavior |
|---|---|---|
| W5 closure | scripts/bqg_v10_w5_lorentzian_adjoint_closure.py | amp_tol defaults to \(10^{-12}\); forward and direct-adjoint columns are filtered before support discovery |
| Ordered KKV accumulation | scripts/peter_weyl_lorentzian_ordered_triple_gate.py | TOL=\(10^{-11}\); sparse additions discard small coefficients |
| Covariant composition | scripts/peter_weyl_covariant_composition_gate.py | TOL=\(10^{-11}\); projection and output maps discard small coefficients |
| Covariant K composition | scripts/peter_weyl_covariant_K_composition_gate.py | TOL=\(10^{-11}\); it calls thresholded composition and complete-state routines |
| Covariant volume leg | scripts/peter_weyl_covariant_volume_leg_gate.py | charged projection uses coefficient cutoffs including \(10^{-13}\) and \(10^{-11}\) |
| Earlier K block | scripts/peter_weyl_lorentzian_K_block_gate.py | local and state pruning thresholds include \(10^{-9}\) and \(10^{-10}\) |
| Volume zero convention | scripts/peter_weyl_zeroaware_volume_migration_experiment.py | eigenvalues below a backward-error threshold are set to zero before \(\sqrt{|Q|}\); this is a separate zero-space convention and must be audited on its own |

The v10 runner's local sparse accumulator uses zero threshold, but that does not remove the thresholds inside the imported Peter-Weyl components. The v10 amp_tol is an additional layer.

The algebraic angular-momentum selection rule that removes source \(J=2\) after the middle leg is different: if the final leg has rank at most one and the target must be a scalar, that forbidden sector can be excluded by an exact Clebsch-Gordan selection rule. It must not be conflated with coefficient-magnitude pruning.

## Exact counterexample

For \(P_\tau(z)=z\) when \(|z|>\tau\), and zero otherwise, take \(\tau=1\), \(x=1.5\), and \(\lambda=0.5\):
\[
P_1(x)=1.5,\qquad P_1(\lambda x)=0,\qquad \lambda P_1(x)=0.75.
\]
Thus the threshold map violates homogeneity.

Accumulation can also depend on grouping. With threshold one and terms \((0.75,0.75,-0.40)\), pruning after every addition gives zero, while summing first and pruning once gives \(1.10\). The standalone regression script reproduces both counterexamples and the same homogeneity defect at a representative \(10^{-11}\) threshold.

## Consequence for the current W5 status

1. A v10 PASS_ADJOINT_CONSISTENCY would establish agreement between the two finite, thresholded constructions on the discovered component. It would not establish exact adjoint closure of the unthresholded KKV map.
2. A cutoff comparison at \(J_{\max 2}=10,12,14\) is not sufficient if the amplitude-pruning error is uncontrolled. Cutoff and amplitude tolerance must be varied independently.
3. The existing random-matrix unit tests validate the reversed-order adjoint index algebra. They do not run the actual W5 Peter-Weyl map and do not test its linearity under rescaling.
4. The historical partial raw-column run (one of 24 epsilon terms completed, then stopped) remains incomplete evidence. It is neither a zero result nor a nonzero result for the full W5 KKV column.

## Required repair before an exact W5 operator claim

### A. Separate structural zeros from numerical smallness

- Keep exact Clebsch-Gordan/triangle/parity selection rules as symbolic support rules.
- Do not remove a coefficient merely because its floating-point magnitude is below a fixed threshold.
- Accumulate all contributions to a matrix element before applying any reporting/display threshold.
- Treat zero-aware spectral truncation of \(Q\) as a separately specified convention, with eigenvalue-gap and backward-error certificates.

### B. Build a linear operator representation

Preferred route: construct matrix elements column-by-column in a fixed orthonormal basis, preserving every coefficient, with a deterministic basis ordering. Cache immutable component matrices keyed by the full discrete labels and exact input basis state; caching may reuse results but must not change coefficients or support.

If exact algebraic coefficients are practical, use exact CG/surd arithmetic. Otherwise use arbitrary-precision ball/interval arithmetic and certify whether each coefficient is zero or provide an explicit unresolved interval. A floating-point threshold alone cannot certify a structural zero.

### C. Validate homogeneity and additivity

For representative basis vectors and superpositions, test
\[
B(\lambda x)=\lambda B(x),\qquad B(x+y)=B(x)+B(y)
\]
at several scales and phases. Report absolute and relative residuals with propagated numerical error bounds. These are necessary numerical checks, not substitutes for an exact construction.

### D. Only then close the adjoint and physical sector

1. Build the complete forward matrix on the structurally generated basis.
2. Independently build the reversed-order adjoint matrix.
3. Check matrix-adjoint residual and leakage.
4. Repeat with increasing \(J_{\max 2}\), precision, and certified coefficient intervals.
5. Select the physical Hermitian ordering and normalization from the operator definition, not from a desired kernel.
6. Construct the non-vacuum master kernel only after these gates pass.

## Reproducible regression

Run: python scripts/bqg_w5_threshold_linearity_counterexample.py

Expected status: PASS_COUNTEREXAMPLES. This test demonstrates the mathematical defect of hard thresholding; it does not execute or certify the expensive W5 physical operator.

## Claim boundary

Established by this audit: the inspected code paths contain coefficient-threshold pruning; hard threshold pruning is not a linear map on the full vector space; therefore thresholded closure cannot be promoted to exact operator closure without a separate error certificate.

Not established: that the underlying formal BQG operator is nonlinear; that the full W5 map is zero or nonzero; that a physical Hermitian Lorentzian constraint exists; or that BQG has a nontrivial physical kernel, HDA, TT sector, or Einstein infrared limit.
