# BCQG IR universality closure theorem

Status: **conditional theorem + exact algebraic reductions; not a replacement for the outstanding direct collective HDA/refinement falsifier.**

## 1. Purpose

The historical collective GR killer lists several simultaneous targets:

\[
D_{space}\to3,\qquad c_{DW}\to\frac12,\qquad z\to1,
\]

\[
(r_G,r_D,r_H,r_{extra})\to(3,3,1,0),\qquad r_{SC}\to0,
\]

\[
N_{phys}\to2,\qquad \Delta_{HH}^{collective}\to0.
\]

These are excellent independent numerical cross-checks, but they are **not logically independent physical miracles**. Under explicit IR hypotheses, several are consequences of the others. This note isolates the genuinely independent burden of proof and prevents redundant brute-force calculations from being mistaken for separate theory axioms.

## 2. Hypotheses

Consider a BCQG refinement family with collective metric carrier and effective scalar constraint. Assume, in one common regulator-safe window:

**H1 — nondegenerate three-dimensional metric phase.** The recursive PL carrier is a 3-manifold and the dynamically reconstructed collective spatial metric remains uniformly nondegenerate. More precisely, on every sufficiently fine level there are constants independent of level

\[
0<c_-\le c_+<\infty
\]

such that the collective metric is locally bilipschitz-equivalent to the intrinsic PL metric.

**H2 — local two-derivative leading IR scalar.** After integrating out gapped microscopic sectors, the leading relevant/marginal scalar constraint is local and of ADM metric form up to higher-derivative/irrelevant operators,

\[
H[N]=\int d^3x\,N\left[
A\frac{\pi_{ab}\pi^{ab}-c\pi^2}{\sqrt q}
-B\sqrt q(R-2\Lambda)
\right]+O(\partial^4/\Lambda_{UV}^2).
\]

**H3 — first-class collective closure.** The independent collective generator ranks stabilize to

\[
(r_G,r_D,r_H,r_{extra})=(3,3,1,0),
\]

with no reducibility pathology in the physical tangent space, and the Hamiltonian commutator approaches the standard hypersurface-deformation target,

\[
[H[N],H[M]]\to i\hbar D[q^{ab}(N\partial_bM-M\partial_bN)].
\]

**H4 — no surviving second-class sector.** The reduced collective constraint matrix has

\[
r_{SC}=0.
\]

The theorem below is conditional on H1--H4. None of these hypotheses may be inserted by target fitting.

## 3. Dimension is not an independent fit once H1 holds

The recursive BCQG carrier is PL three-dimensional. A uniformly nondegenerate bilipschitz change of metric preserves local Hausdorff dimension and the principal second-order elliptic structure of the Laplace operator. Therefore the collective metric phase cannot change its local spatial dimension while H1 remains valid:

\[
\boxed{D_{space}=3.}
\]

A heat-kernel/spectral-dimension calculation remains an excellent independent finite-size falsifier, but it is not a logically separate axiom once uniform metric regularity has been established on the refinement family.

The practical microscopic test is therefore stronger and cleaner than choosing a diffusion-time window: demonstrate level-independent lower and upper singular-value bounds for the flux-to-metric Jacobian / collective metric map.

At the first canonical coarse block the repository already has

\[
M_{hq}^TM_{hq}=\frac14I-\frac1{12}O_{opp},
\]

hence

\[
\boxed{s_{min}=1/\sqrt6,\qquad s_{max}=1/\sqrt3,\qquad \kappa=\sqrt2.}
\]

The open task is to control these bounds under refinement, not to fit dimension independently at every level.

## 4. HDA fixes the DeWitt trace coefficient and relative kinetic/curvature normalization

For the local ADM family in H2 the exact bracket identity already established in the repository is

\[
\{H[N],H[M]\}
=AB\left[D[\beta]+4\left(c-\frac12\right)I[N,M]\right],
\]

where

\[
\beta^a=q^{ab}(N\partial_bM-M\partial_bN).
\]

For generic canonical data and arbitrary lapses, H3 implies

\[
\boxed{c=\frac12,\qquad AB=1.}
\]

Thus a direct measured DeWitt Hessian remains a powerful blind cross-check, but under H2+H3 its target value is a theorem-level consequence of collective HDA closure rather than an additional tunable requirement.

The remaining two-derivative classical freedom is the familiar overall Newton normalization and cosmological constant,

\[
A=16\pi G,\qquad B=(16\pi G)^{-1},
\]

with \(\Lambda\) absent from the HDA bracket.

## 5. Relativistic scaling follows from the same closure

In the transverse-traceless sector,

\[
H_{TT}=\int d^3x\left[A\pi_{TT}^2+\frac B4(\partial h_{TT})^2\right]
\]

gives

\[
\ddot h_{TT}=AB\nabla^2 h_{TT}.
\]

Since H3 forces \(AB=1\),

\[
\boxed{\omega^2=k^2+O(k^4/\Lambda_{UV}^2),\qquad z=1}
\]

in the leading local IR theory. Hence \(z\to1\) is another independent numerical diagnostic, but not an independent low-energy parameter once H2+H3 hold.

## 6. The graviton mode count follows from the constraint ranks

The canonical SU(2)-connection/flux phase-space convention has 18 real phase-space dimensions per generic bulk point/block. With H3 and H4,

\[
r_{FC}=3+3+1=7,
\qquad r_{SC}=0.
\]

Dirac counting gives

\[
N_{phys}^{config}
=\frac{18-2r_{FC}-r_{SC}}2
=\frac{18-14}{2}
=\boxed{2}.
\]

Therefore \(N_{phys}=2\) is a derived check of the rank theorem, not a separate adjustable target. A direct spectrum calculation remains valuable as a held-out falsifier for hidden reducibility or missed second-class directions.

## 7. Canonical depth-two effective scalar: no free denominator

The first homogeneous collective carrier obeys the exact symmetry result

\[
P C P=0,
\qquad P=W_gW_g^\dagger,
\]

for the direct physical Euclidean/Lorentzian gravitational block. Let \(Q=1-P\) after separately retaining any exact low-energy zero modes that must be promoted into \(P\). Write the full constraint on \(P\oplus Q\) as

\[
C=\begin{pmatrix}
0 & B\\
B^\dagger & D
\end{pmatrix}.
\]

On the coupled subspace where \(D=QCQ\) is invertible, the exact constraint equation \(C\Psi=0\) gives

\[
\psi_Q=-D^{-1}B^\dagger\psi_P
\]

and therefore the unique Schur/Feshbach reduction at zero constraint energy,

\[
\boxed{
C_{eff}(0)=-B D^{-1}B^\dagger
=-PCQ\,(QCQ)^{-1}\,QCP.
}
\]

So the depth-two denominator is **not a fit parameter** and should not be selected to improve agreement with GR.

If \(D\) has an exact kernel:

1. components of the kernel coupled to \(P\) must be promoted into an enlarged low-energy carrier;
2. only after exact zero modes are separated may the Moore--Penrose inverse be used on the orthogonal gapped range;
3. no ad-hoc \(i\eta\), mass shift, or denominator clipping may be introduced after consulting the GR target.

This reduces the collective depth-two calculation to the Krylov-reachable coupled sector rather than the full Hilbert space.

## 8. IR universality consequence

Under H1--H4, and within the local two-derivative ansatz H2,

\[
\boxed{
\text{BCQG}_{IR}
\in
\text{ADM/GR}_{G,\Lambda}
+O(\partial^4/\Lambda_{UV}^2).
}
\]

The following quantities are then linked consequences:

\[
\boxed{
D_{space}=3,
\quad c_{DW}=\frac12,
\quad AB=1,
\quad z=1,
\quad N_{phys}=2.
}
\]

This is a universality statement, not an experimental confirmation of BCQG and not a claim that H1--H4 have already all been established dynamically.

## 9. What is now genuinely independent

For internal closure of the gravity sector, the irreducible remaining burden is much smaller than the historical checklist:

1. **metric regularity under refinement:** prove uniform nondegeneracy/bilipschitz control of the dynamically reconstructed collective metric;
2. **collective first-class algebra:** show the frozen collective \([H,H]\) converges to the diffeomorphism target on the same family;
3. **constraint completeness:** establish stable ranks \((3,3,1,0)\) and \(r_{SC}=0\);
4. **IR locality/relevance:** show non-ADM extra operators are gapped or higher-derivative/irrelevant rather than additional unsuppressed low-energy fields;
5. **absolute scale:** match the remaining Newton normalization \(G\) and \(\Lambda\) to microscopic observables.

Everything else in the old collective killer remains valuable as independent regression/cross-validation, but not as a logically independent closure condition.

## 10. Relation to existing repository results

This theorem composes, without changing their scientific status:

- `GLOBAL_MANIFOLD_Q2_COMPLETION.md` and the recursive PL-S3 carrier;
- `COLLECTIVE_METRIC_CALIBRATION_IRREP_THEOREM.md` and its exact first-block singular spectrum;
- `COLLECTIVE_GRAVITATIONAL_DIRECT_BLOCK_THEOREM.md` and \(PCP=0\);
- `DEWITT_HDA_UNIQUENESS.md`;
- `ADM_HDA_PARAMETER_SELECTION.md`;
- `COLLECTIVE_CONSTRAINT_RANK_PROTOCOL.md`;
- `BINARY_HDA_UNIQUENESS_CHAIN.md`;
- the microscopic operator-first HDA architecture certificates.

The preregistered `COLLECTIVE_GR_UNIVERSALITY_PREREGISTRATION.md` remains active as a stronger numerical AND-gate and must **not** be relaxed after seeing data.

## 11. Scientific status

**Proved algebraically under stated hypotheses:** HDA \(\Rightarrow c=1/2,AB=1\); ranks \((3,3,1,0),r_{SC}=0\Rightarrow N_{phys}=2\); exact Schur/Feshbach reduction for the direct-block-zero carrier.

**Conditional geometric consequence:** uniform bilipschitz nondegenerate metric refinement on recursive PL-S3 \(\Rightarrow D_{space}=3\).

**Still open dynamically:** establishment of H1--H4 on the same direct BCQG refinement family, corrected finite Lorentzian collective amplitudes, collective depth-two effective scalar, absolute scale matching, Lorentzian quantum measure/unitarity, complete chiral matter sector, and experiment.
