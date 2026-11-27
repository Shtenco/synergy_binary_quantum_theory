# Cross-branch proof reconciliation: Lorentzian adjoints and causal cones

Date: 2026-11-27

Status: **source-level reconciliation of the repository's branch-local theorems and finite certificates. This note does not promote conditional or finite results into an unconditional physical BQG theorem.**

## 1. Search coverage

The repository currently exposes 167 branch refs (100 on the first page and 67 on the second). Their recursive Git trees were inspected for causal-cone, null-cone, Lorentzian, adjoint, Hermitian, physical-history and metric artifacts. The relevant source documents were then read on the branches that contain them.

Important result: several material proofs are branch-local and are not integrated into the current `main` line. In particular, `audit/bqg-independent-proof-20261108`, `audit/pl48-spectral-complete-20261110`, `audit/v26-strict-leakage-20261110`, and `research/physical-pole-universality` are divergent from `main`. The current `research/null-cone-physical-bridge-v2-20261127` branch is ahead of `main` and contains the latest null-cone audit, but it does not yet integrate the branch-local Hermitian-completion theorem.

## 2. What is proved about adjoints and Hermiticity

### 2.1 The K commutator is anti-Hermitian under its operator hypotheses

For Hermitian (V_v) and Hermitian sine-ordered Euclidean constraint (H^{sine}_{E,v}), on a common invariant domain,

[
K_{v,comm}=[V_v,H^{sine}_{E,v}],
qquad
K_{v,comm}^{\dagger}=-K_{v,comm}.
]

The Hermitian K conversion is

[
K_{v,H}=(i\hbar)^{-1}K_{v,comm}.
]

For two such K legs the factor is exactly (-\hbar^{-2}). This is a statement about the K commutator legs; it is not a statement that the complete unsymmetrized Lorentzian KKV word is anti-Hermitian.

### 2.2 The unsymmetrized raw Lorentzian KKV operator is not globally anti-Hermitian

`LORENTZIAN_HERMITIAN_COMPLETION.md` on `audit/bqg-independent-proof-20261108` records exact fixed-environment blocks with real Pauli (X/Z) components before environment tracing. Those blocks rule out the global identity

[
L_{raw}^{\dagger}=-L_{raw}
]

for the unsymmetrized microscopic ordering.

The physical candidate is instead the Hermitian completion

[
S=-\frac{i}{2}(L_{raw}-L_{raw}^{\dagger}).
]

At (\beta=\hbar=1), the documented signed correction is

[
H_{corr}^{sym}=\frac{16i}{9}(L_{raw}-L_{raw}^{\dagger}).
]

### 2.3 Exact direct-Hermitian commutator theorem exists

`LORENTZIAN_DIRECT_HERMITIAN_COMMUTATOR_THEOREM.md` on `audit/bqg-independent-proof-20261108` and sibling branches proves an exact identity for the frozen v1.3 Lorentzian word. Given the component adjoint identities (C_a(K)_{ij}^{\dagger}=-C_a(K)_{ji}) and (C_c(V)_{ij}^{\dagger}=C_c(V)_{ji}), the six epsilon-ordered words reduce exactly to three commutator words:

[
L_d=\eta_d\sum_{(a,b,c)\in cyclic}
\operatorname{Tr}_{aux}([C_a(K),C_b(K)]C_c(V)).
]

The physical Hermitian block can therefore be written directly as

[
S=-\frac{i}{2}\sum_d\eta_d\sum_{(a,b,c)\in cyclic}
\operatorname{Tr}_{aux}
\left(\{[C_a(K),C_b(K)],C_c(V)\}\right).
]

Each ([C_a(K),C_b(K)]) is anti-Hermitian; its anticommutator with Hermitian (C_c(V)) is anti-Hermitian; multiplication by (-i/2) makes each term Hermitian. This is an exact operator-algebra identity for the declared word and assumptions, without a tetrahedral-covariance assumption.

### 2.4 The finite-cutoff Hermitian projection is uniquely characterized

`LORENTZIAN_HERMITIAN_PROJECTION_UNIQUENESS.md` proves that on a finite Peter–Weyl cutoff block,

[
\operatorname{End}(\mathcal H)=\operatorname{Herm}(\mathcal H)
\oplus\operatorname{AntiHerm}(\mathcal H)
]

as real Hilbert spaces, and that the unique projection onto the Hermitian part is the Hilbert–Schmidt closest Hermitian operator. For (A=-iL_{raw}),

[
P_H(A)=-\frac{i}{2}(L_{raw}-L_{raw}^{\dagger}).
]

The committed finite certificate reports all checks passed, zero closest-point identity defect, and unitary-covariance defect (1.3437330826322668\times10^{-15}). The Hermitian-completion result JSON also reports all its declared checks passed.

### 2.5 Boundary of the operator result

These findings correct any blanket statement that “nothing has been proved about the full Lorentzian operator.” There is a genuine exact algebraic Hermitian-completion theorem and a finite-cutoff uniqueness certificate.

They do **not**, by themselves, establish that the current production Peter–Weyl implementation is already the completed physical constraint on its full reachable domain. The direct-Hermitian theorem's own promotion contract still requires:
- exact symbolic/random-matrix identity against the original forward-plus-adjoint word;
- finite Peter–Weyl comparison against the historical collector;
- a direct 16-cell physical-pair comparison;
- unchanged (V_{tet}), (J_{max,2}=7), zero-aware spectrum, scalar-projection and leakage guards;
- production provenance identifying the direct-Hermitian implementation.

The full-HDA preregistration `PETER_WEYL_HERMITIAN_FULL_HDA_PREREGISTRATION.md` also says a completed finite HDA run must report Hermiticity of all represented (S_0,S_1,G_0,G_1) blocks and the signed channel Gram. It explicitly distinguishes a finite PASS from a continuum or experimental proof.

The W5 adjoint-closure audit is a separate result: it specifies how to compare the raw forward KKV matrix with its independently ordered adjoint on a completely closed finite reachable component. Its own document says that even a PASS would not select the physical Hermitian ordering, restore the physical prefactor, establish cutoff convergence, or prove a physical joint kernel. The W5 production pipeline audit still records those steps as blockers.

## 3. What is proved about causal/light cones

The branch inventory contains multiple exact cone-related results. They are different mathematical objects.

### 3.1 Exact reduced-model propagation cone

`BQG_SHARP_ANALYTIC_CAUSAL_CONE.md` proves the sharp exponential propagation front for the homogeneous reduced free-fermion chain and its gapped staggered-volume deformation:

[
v_{front}(m)=\sqrt{4m^2+9}-2|m|.
]

This is an exact theorem within that reduced model. The same document explicitly says it is not the physical speed of light of full BQG and does not establish Lorentz invariance, universality across species, or a continuum metric null cone.

### 3.2 Exact graph-sector causal-order bounds

`BQG_GRAPH_CHANGING_CAUSAL_ORDER_THEOREM.md` proves finite-regulator graph-sector support/order bounds for powers of the master constraint and local relational steps. It also explicitly leaves open the identification of master-flow parameter (T) with physical spacetime time, and the full Lorentzian spacetime causal cone.

### 3.3 Conditional ADM/HDA tensor cone

The ADM/HDA parameter-selection result gives (c=1/2), (AB=1), and (c_T=1) inside the declared local two-derivative ADM family. This is a conditional leading-TT characteristic target. The separate physical TT kernel document states that the theory-specific connected physical (\Gamma^{(2)}_{TT,phys}) remains open.

### 3.4 Common photon–graviton light-cone theorem exists, but is conditional

`COMMON_PHOTON_GRAVITON_LIGHT_CONE_THEOREM.md` on `research/physical-pole-universality` proves that the leading photon and graviton cones coincide if the candidate theory reaches an infrared phase with:
1. one emergent Lorentzian metric;
2. an unbroken massless metric spin-2 sector;
3. a deconfined compact (U(1)) sector with leading Maxwell action;
4. both sectors coupled to the same metric/Hodge geometry;
5. a locally Lorentz-invariant vacuum with no surviving preferred tensor or foliation order parameter.

Under these hypotheses both principal symbols have the same metric null equation. This is a valid conditional infrared effective-theory theorem. It is not, by itself, a derivation from the microscopic BQG operators that all five hypotheses are realized.

The same branch's `PHYSICAL_CLOSURE_LEDGER_V3.md` labels gravity/light-cone universality as a conditional IR theorem and explicitly lists the theory-specific physical-history kernel, connected physical metric correlators, physical graviton pole and photon stiffness/deconfinement as open physical gates.

### 3.5 Latest null-cone bridge CI is a regression pass, not a physical-cone derivation

On `research/null-cone-physical-bridge-v2-20261127`, workflow run `46040510729` passed 1,521 algebraic regression checks:
- 500 ADM null-coordinate identity checks;
- 9 HDA/TT speed-normalization checks;
- 11 reduced spectral-front checks;
- 1,001 graph support-order checks.

The run log itself ends with: `SCOPE: algebraic bridge regression only; no physical BQG null cone derived`. The gate script has the same explicit scope statement.

The latest `BQG_PHYSICAL_NULL_CONE_BRIDGE_AND_TWO_INPUT_AUDIT_20261127.md` therefore correctly says that the full physical BQG Lorentz cone has not yet been derived by the currently integrated operator chain.

## 4. Correct reconciled statement

The strongest defensible statement is:

1. **K-leg anti-Hermiticity:** exact under the Hermiticity/common-domain hypotheses for (V_v) and (H^{sine}_{E,v}).
2. **Full Lorentzian raw word:** not globally anti-Hermitian; fixed-environment finite blocks contain a nonzero Hermitian part.
3. **Physical Lorentzian completion:** an exact direct-Hermitian operator identity exists for the frozen v1.3 word; the finite-cutoff minimal Hermitian projection is uniquely characterized and its saved certificate passes.
4. **Reduced and graph causal cones:** exact within their declared finite/reduced models.
5. **Common photon–graviton cone:** proved conditionally in the shared-metric Lorentz-invariant infrared EFT.
6. **Unconditional physical null cone of full BQG:** not yet established by the current integrated microscopic physical-history chain. The physical clock/rigging map, connected physical TT kernel, retarded response, common-domain HDA/normal evolution, and controlled refinement/universality gates remain required.

This is not a denial of the branch-local theorems. It is the distinction between an exact theorem under stated hypotheses, a finite certificate, a conditional infrared consequence, and a microscopic physical prediction.

## 5. Required integration work

1. Port the direct-Hermitian identity into the current branch only after matching its frozen v1.3 component definitions to the current sine-K/W5 operator definitions.
2. Verify component adjoints on the same common regulated Peter–Weyl domain, including the complete direct-adjoint ordering and all amplitude-loss layers.
3. Compare the direct-Hermitian producer with the historical forward-plus-adjoint production columns on the declared 16-cell and K5 habitats.
4. Complete the W5 scalar-covariance and full reachable-space adjoint/closure audits without projecting away rejected sectors.
5. Use the completed Hermitian Lorentzian constraint in the common-domain HDA and physical-history/rigging-map construction.
6. Derive the physical relational time or boundary-history ordering and the connected retarded response.
7. Reconstruct (q_{ab},K_{ab},N,N^a) on the same physical states, obtain (g_{\mu\nu}), and compute its null cone.
8. Compare scalar, vector, TT and matter fronts on a common refinement trajectory, independently of the spectral metric calibration.

No parameter fit to the reduced (v_{front}(m)) may substitute for steps 5–8.
