# BQG normal-evolution and extrinsic-curvature operator audit — 2026-11-27

Status: **operator-provenance audit; no microscopic Lorentzian normal generator or physical extrinsic-curvature operator identified in the inspected chain.**

This audit is intentionally conservative. It distinguishes a graph cochain used to label the HDA structure-function support from a lapse/shift field that generates physical normal/tangential deformations. It also distinguishes a finite Euclidean/canonical constraint kernel from Lorentzian relational evolution.

## 1. Executive result

The current inspected implementation contains:

1. A finite K5 Peter–Weyl/Gauss Hilbert space and local Thiemann-type constraint kernels.
2. A graph-changing constraint action and a compressed finite master-constraint shell.
3. A dual-K5 lapse cochain with exact antisymmetric support.
4. A classical local ADM ansatz whose HDA bracket selects (c=1/2), (AB=1).
5. A finite relational-history positive control using an externally specified clock and (R=J).

It does **not** yet provide, in one common physical construction:

- a microscopic lapse operator (\hat N(x)) with a demonstrated continuum/geometric interpretation;
- a shift/vector-deformation operator (\hat N^a(x)) derived from the flux/metric data;
- a self-adjoint Lorentzian normal-evolution generator with a specified domain;
- a physical relational clock derived from the constrained BQG degrees of freedom;
- an operator equation yielding (\dot q_{ab}\) on physical states;
- a derived (\hat K_{ab}\), or a refinement-stable expectation/correlation function for it.

Consequently, the physical Lorentzian null cone is not yet derived by this operator chain. This is a precise missing-link result, not a claim that no such construction can exist.

## 2. Operator dependency and provenance table

| Object | Existing construction | What it establishes | What it does not establish |
|---|---|---|---|
| Link holonomy (U_e) | Finite Peter–Weyl link operators in the K5 implementation | Regulated graph/representation changes | Spacetime time evolution |
| Node volume (V_v) | Active-link projector in the minimal (j=1/2) Thiemann-type gate; other geometric volume constructions are separately regulated | Finite nondegeneracy/volume selection in that gate | A continuum volume density or Lorentzian clock |
| Local constraint (C_v) | (\mathrm{Tr}[(U_{vabv}-U_{vbav})U_{vc}[U_{cv},V_v]]), with finite Gauss projection | Euclidean/canonical finite constraint-kernel tests | A complete Lorentzian Hamiltonian constraint with proven HDA closure |
| Master constraint (\mathbb M=\sum_v C_v^\dagger C_v) | Positive finite shell construction | Kernel and compressed projector for that chosen shell | Physical time; (T) in (e^{-T\mathbb M}) is not automatically time |
| Graph signature | Delete (j=0) links | Abstract graph-sector labels | A spacetime metric or causal cone |
| Lapse cochain (\omega_{vw}=N_vM_w-N_wM_v) | Exact dual-K5 edge identity | Antisymmetric support of a discrete HDA-like structure-function input | A metric-raised shift, physical lapse field, or normal generator |
| ADM (q_{ab},\pi^{ab}) | Classical local ansatz | Conditional closure identities (c=1/2, AB=1) | Microscopic convergence of BQG operators to ADM variables |
| Relational clock (S), system step (R=J) | Exact finite (C_8\) positive control | Mathematical mechanism of clock/system character matching | A clock or evolution operator derived from the BQG Hamiltonian |
| (K_{ab}) | No microscopic operator found in the inspected files | Only the classical target definition below | No BQG-derived extrinsic curvature yet |

## 3. The lapse-cochain result: exact but narrower than a physical shift

For node values (N_v,M_v), define midpoint averages and edge differences on an oriented edge ((v,w)):
[
\bar N_{vw}=\frac{N_v+N_w}{2},\qquad
(dM)_{vw}=M_w-M_v,
]
and similarly for (M). Then
[
\omega_{vw}=\bar N_{vw}(dM)_{vw}
-\bar M_{vw}(dN)_{vw}
=N_vM_w-N_wM_v.
]
The identity is exact. For basis lapses (N=\delta_i), (M=\delta_j), only the shared edge ((i,j)) is supported.

This is the correct combinatorial antisymmetric cochain associated with the discrete lapse pair. It is **not yet**
[
\beta^a=q^{ab}(N\partial_bM-M\partial_bN),
]
because the geometry-dependent sharp map from the dual edge cochain to a spatial vector field has not been supplied. In particular, it cannot be called the physical shift operator without:
1. a nondegenerate metric/flux map on the same state space;
2. a demonstrated discrete differential and inverse/sharp map with controlled kernel;
3. covariance under the relevant graph/gauge transformations;
4. a refinement limit reproducing the continuum structure function.

The lapse cochain is therefore a structural input to the HDA test, not the output of that test.

## 4. Extrinsic curvature: required derivation and sign convention

Fix signature ((-+++)), future-directed unit normal (n^\mu n_\mu=-1), and
[
K_{ab}:=-q_a{}^\mu q_b{}^\nu\nabla_\mu n_\nu.
]
With
[
t^\mu=Nn^\mu+N^\mu,
]
this convention gives
[
\boxed{K_{ab}=\frac{1}{2N}
(\dot q_{ab}-\mathcal L_{\vec N}q_{ab}).}
]
The opposite convention for (K_{ab}) changes the overall sign and must not be mixed with formulas from another source.

For BQG, the right-hand side is not an operational definition until the following are independently specified:

- (\hat q_{ab}) (or a coarse-grained relational metric observable) on a common physical domain;
- the physical relational parameter (\tau) and clock observable (T_{\rm clock});
- the generator (\hat H_{\rm rel}) that evolves relational observables, including its self-adjointness/domain or a well-defined boundary-amplitude replacement;
- the tangential deformation (\hat{\mathcal L}_{\vec N}\hat q_{ab}), with (\vec N) derived from the same geometry;
- the lapse normalization, so that the normal is unit normalized in the same emergent metric.

Only then can one define, for example,
[
\frac{d\hat q_{ab}}{d\tau}
= i[\hat H_{\rm rel},\hat q_{ab}]
]
(in units (\hbar=1\), when a genuine relational Hamiltonian exists), and test the classical-limit identity for (K_{ab}). A projector (P_{\rm phys}) alone is insufficient: it selects a constraint sector but does not specify a unique relational flow.

## 5. Why current candidates do not close this gap

### 5.1 K5 Thiemann-type gate

The current K5 kernel is a finite Euclidean/canonical construction. Its own scope note explicitly withholds the full real Lorentzian Hamiltonian, HDA closure, continuum universality, and matter/experiment. Its common-kernel result is valuable but is not a normal-evolution generator.

### 5.2 Oriented quantum-HDA diagnostic

The oriented (J_{\max}=1/2) diagnostic restores permutation covariance but explicitly fails the physical HDA gate: the cutoff is below the Peter–Weyl HH safe wall and the fixed-sector commutator has large graph leakage. It cannot be promoted to a physical shift or (K_{ab}).

### 5.3 Compressed graph-changing master shell

The shell uses actual regulated (H_v) actions, but compresses them to a finite generated basis. Its own scope requires depth/cutoff convergence and a stable zero-mode projector before interpreting cross-graph blocks. Moreover, the master parameter (T) is constraint/projector flow, not physical time.

### 5.4 Relational (C_8) positive control

The clock shift and (R=J) are externally declared positive-control objects. The construction proves the algebraic Page–Wootters mechanism, while explicitly not deriving a physical clock or microscopic gravitational evolution operator.

### 5.5 Classical ADM closure

The (c=1/2, AB=1) result is conditional on the local ADM ansatz. It tells us what a successful BQG continuum limit must reproduce; it does not supply the missing microscopic normal generator.

## 6. Required closure sequence (no parameter fitting)

The next physical program should be ordered as follows:

1. **Microscopic Lorentzian constraint:** specify the actual Lorentzian (\hat C[N]), including its regulator, adjoint, domain, and graph-changing action. Do not infer it from the Euclidean K5 kernel.
2. **Constraint algebra:** evaluate the regulated commutators on a common reachable Peter–Weyl habitat and quantify graph leakage, Gauss projection leakage, and regulator dependence.
3. **Structure-function map:** construct the flux/metric-dependent discrete sharp map taking the lapse cochain to the vector deformation; prove covariance and convergence to (q^{ab}(N\partial_bM-M\partial_bN)).
4. **Relational clock:** derive a clock observable from the same microscopic theory, or give a boundary-amplitude construction with a controlled relational parameter. Do not insert (R=J) as physical evolution.
5. **Metric evolution:** compute relational changes/correlators of (q_{ab}) using the physical generator/amplitude.
6. **Extrinsic curvature:** construct (K_{ab}) with the fixed sign convention and test its symmetry, gauge covariance, classical limit, and refinement stability.
7. **Lorentzian metric:** assemble (g_{\mu\nu}) from (q_{ab},N,N^a), with unit-normal normalization fixed.
8. **Null cone:** only now solve (g_{\mu\nu}v^\mu v^\nu=0), then feed the resulting null generators into the already implemented Jacobi/Riccati optical gate.

A failure at steps 1–4 blocks any claim that the later cone is a physical BQG prediction.

## 7. Minimal acceptance tests

A candidate closure is not accepted unless all of these pass in one common regulator/refinement window:

- **Domain:** all operators act on a stated common dense/finite regulated domain; adjoints and boundary terms are controlled.
- **Gauss/diffeomorphism covariance:** constraint and metric observables preserve the required gauge sectors.
- **HDA:** ( [\hat H[N],\hat H[M]] ) approaches the diffeomorphism generator with the same geometry-derived structure function as the independent lapse-cochain/sharp-map construction.
- **No anomaly by projection artefact:** report both full reachable-space commutator and projected-sector result; quantify discarded norm.
- **Clock independence:** physical predictions are stable under admissible clock reparameterizations, not tuned to one arbitrary clock step.
- **Extrinsic-curvature identity:** the commutator/relational derivative of (q_{ab}) reproduces the stated (K_{ab}) relation, including shift subtraction and lapse normalization.
- **Classical limit:** (c_{eff}\to1/2), (A_{eff}B_{eff}\to1), and the TT characteristic speed tends to one in normal units without retuning.
- **Refinement:** all preceding residuals converge as graph/spin cutoffs are increased; no finite-shell zero mode is silently treated as a continuum physical state.

## 8. Current conclusion

The exact dual-K5 lapse cochain is a real piece of the discrete HDA structure, and the classical ADM/HDA chain fixes a sharp target. But the inspected microscopic chain still stops before a physical Lorentzian normal generator. The next valid milestone is **not** to draw or fit another cone. It is to derive the Lorentzian constraint plus its relational generator, then construct (K_{ab}) from the same (q_{ab}), lapse, and shift.

Until that milestone is met, the optical Jacobi/Riccati layer remains a valid downstream mathematical gate, but its input null congruence is not yet established as a physical BQG null congruence.
