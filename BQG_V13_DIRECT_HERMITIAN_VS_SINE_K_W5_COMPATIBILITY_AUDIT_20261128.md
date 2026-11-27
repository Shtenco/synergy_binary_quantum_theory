# v1.3 direct-Hermitian theorem vs sine-K/W5 production: compatibility audit

Date: 2026-11-28

Status: **operator-level compatibility is partial; production promotion is blocked by unresolved phase bookkeeping and unclosed finite-domain adjoint tests.** This audit does not alter the production operator, thresholds, normalization, or main branch.

## 1. Executive finding

The v1.3 theorem and current W5 code share the same *candidate K-commutator algebraic pattern*:

- v1.3 assumes each matrix-covariant K component obeys
  `C_a(K)_{ij}^dagger = -C_a(K)_{ji}`;
- current sine-K constructs
  `H_E^sine=(T-T^dagger)/(2i)`, then `K_comm=[V,H_E^sine]`;
- if V and H_E^sine are Hermitian on one common invariant domain, then K_comm is anti-Hermitian.

That is not yet enough to port the theorem into production. Two independent issues remain:

1. the current finite W5 implementation has not established the component adjoint identities on its complete reachable basis with all nested cutoffs accounted for;
2. the W5 normalization ledger currently defines a *different Hermitianization* from the v1.3 direct-Hermitian theorem.

The second issue is algebraic, not merely a missing numerical test. It must be resolved before any implementation is selected as the physical Lorentzian operator.

## 2. Exact v1.3 statement

The frozen v1.3 word is

[
L_{raw}=sum_d\eta_d\sum_{(a,b,c)\in cyclic}
\operatorname{Tr}_{aux}([C_a(K),C_b(K)]C_c(V)).
]

Under the component adjoint assumptions K anti-Hermitian and V Hermitian, the theorem gives

[
S_{v1.3}=-\frac{i}{2}(L_{raw}-L_{raw}^\dagger)
=-\frac{i}{2}\sum_d\eta_d\sum_{cyclic}
\operatorname{Tr}_{aux}\{[C_a(K),C_b(K)],C_c(V)\}.
]

The proof does not use cyclicity of the partial auxiliary trace. It uses the exact reversal of operator order under adjoint and the two component-adjoint identities. The output is Hermitian term by term.

Source: `LORENTZIAN_DIRECT_HERMITIAN_COMMUTATOR_THEOREM.md` on `audit/bqg-independent-proof-20261108`.

## 3. Current sine-K/W5 definitions actually present in code

### 3.1 Euclidean and K operators

[
H_E^{sine}=\frac{T-T^\dagger}{2i},\qquad
K_{comm}^{sine}=[V,H_E^{sine}].
]

The adapter `scripts/peter_weyl_covariant_K_sine_composition_gate.py` installs this sine-Hermitian operator in both the complete-charge and Gauss representations, records `old_plus_K_reused=False`, and applies the historical Gauss K cutoff (10^{-9}).

The W5 normalization ledger correctly records the conditional operator fact

[
K_{comm}^\dagger=-K_{comm}
]

when V and H_E^sine are Hermitian on a common invariant domain. It also records the Hermitian-K conversion

[
K_H=(i\hbar)^{-1}K_{comm},\qquad
K_{H,a}K_{H,b}V_c=-\hbar^{-2}K_{a,comm}K_{b,comm}V_c.
]

Source: `BQG_W5_K_COMMUTATOR_NORMALIZATION.md`.

### 3.2 Matrix-covariant legs and raw W5 map

The current W5 path uses `C_e(V)` and `C_e(K_sine)`, then computes the ordered triple

[
T_{abc}=\operatorname{Tr}_{aux}[C_a(K_{comm})C_b(K_{comm})C_c(V)]
]

with right-to-left state action. The v10 closure runner independently builds the reversed-order adjoint

[
T_{abc}^\dagger=\operatorname{Tr}_{aux}[C_c(V)C_b(K_{comm})C_a(K_{comm})],
]

including transposed auxiliary indices. It explicitly avoids partial-trace cyclicity.

The raw node map is the coherent 24-term epsilon sum over the four omitted-edge faces and six permutations per face. This is structurally the same K-K-V word family as the v1.3 theorem, provided the exact C-leg definitions, edge orientation, epsilon signs, and auxiliary-index convention are matched.

Sources:
- `scripts/peter_weyl_lorentzian_sine_ordered_triple_gate.py`
- `scripts/bqg_v10_w5_lorentzian_adjoint_closure.py`
- `BQG_W5_LORENTZIAN_ADJOINT_CLOSURE.md`
- `BQG_LORENTZIAN_EPSILON_KKV_ADJOINT_AUDIT.md`

## 4. Critical normalization / five-bracket phase accounting

The current W5 normalization note gives, after converting both K_comm legs to Hermitian K_H,

[
B_H=-\hbar^{-2}B_{comm},
]

and then defines the symmetric completion

[
B_{H,sym}=-\frac{1}{2\hbar^2}
(B_{comm}+B_{comm}^\dagger).
	ag{W5-plus}
]

The v1.3 direct-Hermitian theorem, however, defines

[
S_{v1.3}=-\frac{i}{2}
(B_{comm}-B_{comm}^\dagger).
	ag{v1.3-minus}
]

These are **not the same operator if the extra five-bracket phase is omitted**. They select different Hermitian components of the same complex structural map.

The exact conversion bookkeeping is:

\[
B_H=-\hbar^{-2}B_{comm},
\]

and hence

\[
\operatorname{Herm}(B_H)
=-\frac{1}{2\hbar^2}(B_{comm}+B_{comm}^\dagger),
\]

whereas restoring an explicit (i) in the five-bracket/Lorentzian prefactor gives

\[
\operatorname{Herm}(iB_H)
=\frac{i}{2}(B_H-B_H^\dagger)
=\hbar^{-2}S_{v1.3}.
\]

Thus the v1.3 direct-Hermitian theorem and the W5 Hermitian-K conversion **are compatible if and only if the full physical prefactor carries the corresponding explicit phase (i)** (with the remaining real normalization handled separately). If the W5 structural map is instead declared to enter with no such phase, the plus-completion is a different operator. The current W5 normalization note separates the additional Lorentzian prefactor from (B_H), but does not yet freeze its phase in the same ledger. That is the concrete bookkeeping gap.

Do not silently replace one formula by the other. Freeze the full regulated classical-to-quantum convention—including the placement of (i), Poisson-bracket factors, beta, kappa, hbar, orientation and trace normalization—before choosing the physical W5 completion. A finite Hermiticity identity alone does not determine that prefactor.

Required decision record:
1. state whether the v1.3 K symbol denotes the anti-Hermitian commutator K_comm or the Hermitian canonical K_H;
2. write the complete five-bracket/Lorentzian prefactor before and after K conversion;
3. derive the resulting operator from the frozen classical constraint, without fitting coefficients to a kernel or HDA residual;
4. only then freeze the physical W5 candidate.

## 5. Compatibility matrix

| Layer | v1.3 theorem | Current sine-K/W5 | Verdict |
|---|---|---|---|
| Euclidean input | frozen v1.3 E/V definitions | (H_E^{sine}=(T-T^\dagger)/(2i)) | candidate match; exact source-definition match still required |
| K leg | anti-Hermitian K commutator | (K_{comm}=[V,H_E^{sine}]) | algebraically compatible under common-domain hypotheses |
| C(K) adjoint | assumed component identity with auxiliary-index transpose | built through charged Peter–Weyl composition and finite cutoff | not yet certified on complete W5 reachable domain |
| C(V) adjoint | Hermitian component identity | covariant volume leg with charged blocks and holonomy closure | must be certified with identical domain and index convention |
| KKV ordering | epsilon K-K-V word | 24-term epsilon K-K-V map | same word family; exact orientation/index map still needs a direct identity test |
| Adjoint | reversed product, no partial-trace cyclicity | independently reversed v10 construction | structurally compatible; numerical closure run still required |
| Hermitian completion | (-i(B-B^\dagger)/2) | W5 note currently specifies (-(B+B^\dagger)/(2\hbar^2)) after K conversion | **incompatible until convention is resolved** |
| Regulator | frozen v1.3 finite cutoff/certificate | W5 seeds (j=2), initial (J_{max,2}=10), nested cutoffs | different finite domains; no certificate transfer |
| Physical prefactor | v1.3 frozen word convention | beta/kappa/hbar/trace factors still separately open in W5 ledger | not matched |
| Continuum/physical cone | not implied by finite theorem | physical history, rigging, retarded response and refinement remain open | not established by this operator audit |

## 6. CI and finite numerical evidence rechecked

### 6.1 Sine-K seed/reverse-row audit

Workflow run `46050985436` (commit `8932377f9068387084e1af8def5358b77f548712`) is **cancelled**. The Python process was interrupted during the seed-column/reverse-row calculation. The log says no `BQG_SINE_K_SEED_ROW_AUDIT.json` was produced; the upload step explicitly reports no file. Thus this run supplies no numerical anti-Hermiticity verdict.

The runner's own docstring correctly limits its claim to the complete reverse row associated with one seed column, not every pair in the full reachable-space matrix. Its call graph also retains internal coefficient/output cutoffs.

### 6.2 W5 scalar tolerance ablation

Workflow run `46051047906` completed successfully as a computation and uploaded artifact `15225270977`. The log reports:
- full-map L2 agreement with the unpruned-`fsum` reference at approximately (1.59\times10^{-15});
- forbidden-sector amplitude at roundoff scale for that diagnostic;
- maximum basis leakage (=1.0);
- support symmetric difference (=898) against the production baseline.

The script/report explicitly states that outer zero tolerances do not disable nested primitive, charged-basis, volume, and state cutoffs. This is useful evidence about outer accumulation sensitivity; it is not a full no-pruning certificate or a W5 scalar-covariance PASS.

## 7. Exact implementation gaps to close

### Gate A — freeze convention before coding

Create a one-page signed convention ledger resolving the mismatch in §4. It must define K_comm, K_H, the phase multiplying the KKV word, the Lorentzian coefficient, and the exact physical completion. No numerical output may be used to choose signs or coefficients.

### Gate B — component adjoints on one shared finite domain

On the same W5-seeded complete reachable basis and at each fixed cutoff:
1. build every (C_a(K_{comm})_{ij}) and compare its independently assembled adjoint with (-C_a(K_{comm})_{ji});
2. build every (C_c(V)_{ij}) and compare its adjoint with (C_c(V)_{ji});
3. preserve all rejected charge/spin sectors and record their squared norm before any scalar projection;
4. record discarded amplitude/norm at every nested threshold (basis overlap, primitive branch, volume block, T-state, map accumulation, closure amplitude);
5. require complete forward-and-adjoint closure, not only one seed reverse row.

### Gate C — prove the exact 12-word identity against the actual 24+24 implementation

For random finite matrices satisfying the declared component adjoint identities, verify

[
-\frac{i}{2}(B-B^\dagger)
=-\frac{i}{2}\sum_{d,cyclic}\eta_d
\operatorname{Tr}_{aux}\{[C_a(K),C_b(K)],C_c(V)\}.
]

Then run the same equality on actual W5 Peter–Weyl component matrices. Compare the direct 12-word output against the independently assembled 24-forward + 24-reversed-adjoint result, entry by entry, before scalar projection.

### Gate D — preserve existing physical guards

Run the exact direct-Hermitian and independently assembled reference through the same:
- W5 seed set (K\in\{0,2,4,6,8\});
- declared (J_{max,2}) walls and next refinements;
- zero-aware spectrum convention;
- scalar projection and full leakage ledger;
- all-node coherent sum and node-lapse convention;
- historical collector comparison.

Do not relax a threshold or discard a support sector to obtain agreement.

### Gate E — only after A–D, integrate HDA and physical chain

Use the convention-frozen completed operator in the full HDA/Gram tests. Then continue to common-domain normal evolution, physical clock/history ordering, rigging map, connected retarded response, reconstruction of (q_{ab},K_{ab},N,N^a), and the physical metric/null cone along a controlled refinement trajectory.

A finite W5 operator identity is a necessary bridge, not the physical spacetime-cone theorem itself.

## 8. Promotion decision

**Do not port the v1.3 formula directly into the W5 production path yet.** The algebraic KKV pattern is promising and the independently reversed v10 adjoint is the correct structural comparison. But the W5 plus-completion and v1.3 minus-completion currently disagree, and the actual finite component-adjoint hypotheses are not yet closed with a complete leakage ledger.

The next valid deliverable is not a new claim of Hermiticity. It is:
1. a frozen convention decision;
2. a component-adjoint certificate on the actual W5 domain;
3. a machine-checked equality between the direct 12-word formula and the independent 24+24 implementation;
4. a same-domain W5 result with no hidden projection or support loss.

Until all four exist, the correct status is **COMPATIBLE K-LEG ALGEBRA / UNRESOLVED PHYSICAL COMPLETION / W5 DOMAIN AUDIT OPEN**.


## 9. Executed convention-discriminator CI

An audit-only deterministic matrix guard was added:

- `scripts/bqg_v13_w5_completion_convention_guard.py`
- `.github/workflows/bqg-v13-w5-completion-convention-guard.yml`

Run `46060510466` completed with `success` and uploaded artifact `15228821370`.

The operator-valued (2\times2) auxiliary / 5-dimensional graph-space witness reported:

| Check | Result |
|---|---:|
| K legs anti-Hermitian | PASS |
| V legs Hermitian | PASS |
| direct 12-word formula vs v1.3 24-forward + 24-adjoint completion | PASS |
| v1.3 completion Hermitian | PASS |
| W5 plus-completion Hermitian | PASS |
| the two prescriptions are distinct on the witness | PASS |

Numerical diagnostics:
- direct 12-word vs v1.3 completion relative error: (3.7520864821232615\times10^{-16});
- v1.3 Hermiticity defect: (0);
- W5 plus-completion Hermiticity defect: (0);
- relative difference between plus and minus prescriptions: (1.5436847137229295).

The first run exposed a JSON serialization bug in the added guard; it was corrected. The final run `46060714237` completed successfully and uploaded artifact `15228516706`.

The final **generic algebraic discriminator** confirms:
- direct 12-word vs v1.3 24-forward + 24-adjoint relative error: (3.6377005628978436\\times10^{-16});
- v1.3 completion Hermiticity defect: (0);
- no-extra-phase W5 plus-completion Hermiticity defect: (0);
- \\(\\operatorname{Herm}(iB_H)\\) equals the v1.3 result;
- no-extra-phase plus vs phase-completed minus relative difference: (1.5436847137229295).

This resolves the *algebraic bookkeeping question*: the v1.3 theorem is compatible with the W5 K_comm-to-K_H conversion when the explicit five-bracket phase (i) is retained. It is not compatible with silently omitting that phase. The guard still does not test actual W5 Peter–Weyl component matrices and does not select the physical prefactor from the classical constraint.

Final run: https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/46060714237

The audit-only script and this document are on `research/null-cone-physical-bridge-v2-20261127`; `main` was not modified.
