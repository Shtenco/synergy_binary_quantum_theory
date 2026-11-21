# Текущий статус BQG — 12 октября 2026

**structural_candidate_closed=false; physicalization_complete=false; experimentally_confirmed=false.** Девять обязательных физических условий остаются открытыми. Исторические декларации `true` ниже сохраняются только как архив и не являются актуальным заключением.

Главный обзор: [README](README.md). [Пути доказательства](PROOF_ROADMAP_RU.md). [Текущий аудит](AUDIT_CURRENT.md). [Инфотон, обе размерностные гипотезы и история](docs/information_infoton_review/REVIEW_RU.md). [Компактные голономии](docs/compact_holonomies/INTEGRATION_RU.md).

Синусный свободный benchmark завершён в своей области: положительный вакуумный спектр и сглаживание с ведущей поправкой a²/(40L²). Из полной физической динамики BQG этот коррелятор не выведен. Полные H_L, согласованный ненулевой предел и interacting TT остаются открытыми. Первоначальный K5 `passed=false` объяснён неустойчивым абсолютным порогом ранга; исправленный нормированный SVD-контроль имеет `passed=true`, ранг 1 и остаток 1/91 в floating point. Это не точный сертификат и не физическая завершённость. [Последнее продолжение](docs/numerical_continuation/README_RU.md) содержит выбранный слабый обратный элемент и начало второй CK-стадии. Модельный инфотонный резонанс не является обнаружением частицы. Переходы в дополнительные пространственные измерения и 2T не доказаны.

## Исторический материал

Ниже сохранены прежние записи. Даты, статусы и формулировки относятся к своим историческим срезам; при противоречии приоритет имеет текущий аудит и исходные свидетельства.

# Theory status — structural candidate versus physicalization frontier

> ## Canonical update — 21 November 2026
>
> The fail-closed finite depth-6 ledger is authoritative: `depth6_frontier.json`, verified by `scripts/verify_depth6_frontier.py`.
>
> Current finite depth-6 status:
>
> - corrected K5 shell: 264,962 Gauss-admissible spin assignments;
> - Hilbert dimension: 3,111,637;
> - four sectors have finite numerical closure certificates: `[1^5]`, `[5]`, `[4,1]`, and `[2,1,1,1]`;
> - the `[2,1,1,1]` sector is CLOSED_FINITE_NUMERICAL: its 16 exact H0-null directions are all lifted by H1 (rank 16/16, sigma_min approximately 1.13045219), and its certified master-kernel dimension is zero;
> - structural support is closed for all seven S5 irreps;
> - three independent numerical master witnesses remain: `[3,2]` (recovery-first; 130,903-dimensional sector), `[3,1,1]` (153,455), and `[2,2,1]` (130,503);
> - the full finite depth-6 theorem `ker M^(d=6)=span{|0>}` is still NOT_YET_PROVED.
>
> The former planned 130007 x 130007 H0-only minor is explicitly superseded; it is not a current blocker. The 16 exact H0-null directions are real, but the preserved H0-complement certificate plus fresh H1 lift closes the `[2,1,1,1]` master sector.
>
> Even completion of all remaining finite sectors will not close the continuum theory. Refinement stability, the theory-specific physical projector/rigging-map limit, connected physical history, the interacting TT kernel, the physical six-Wilson vector and one common physical scale remain open.


**Canonical status update: 2026-09-29**

The repository now records two different scientific propositions separately.

```text
STRUCTURAL / INTERNAL CANDIDATE
structural_candidate_closed              = true
candidate_mathematical_framework         = true

PHYSICALIZATION / OBSERVATION
physical_projector_history_closed        = false
connected_interblock_history_closed      = false
physical_TT_kernel_frozen                = false
IR_six_wilson_vector_frozen              = false
common_physical_scale_calibrated         = false
experimentally_confirmed                 = false
```

These statements are deliberately not collapsed into one `theory closed` Boolean.

`structural_candidate_closed=true` means that the **declared structural candidate construction**, in its registered exact / finite-tested / explicitly conditional scopes, has an end-to-end internal evidence chain.  It does not mean that the theory-specific physical Hilbert space, continuum rigging map, interacting graviton 1PI kernel or experimental truth has been derived.

Machine-readable sources are:

- `theory_gates.json` — structural candidate ledger;
- `physicalization_gates.json` — fail-closed physicalization ledger;
- `CANONICAL_THEORY_PACKAGE.md` — structural evidence index;
- `OPEN_PROBLEMS.md` — open physicalization, stronger extensions and experiments.

---

## 1. Structural chain currently reproduced

The registered structural candidate is

```text
q=2 binary microstructure
 -> exact refinement fixed point d*=3
 -> regular-tetrahedral Walsh flux carrier
 -> face-qubit / Gauss-singlet geometry carrier
 -> graph-changing active/no-link Peter-Weyl representation
 -> selected recursive PL 3-manifold + exact carrier gluing
 -> B / simplicity / Urbantke metric / compatible connection
 -> Regge / Einstein-Hilbert controls
 -> ADM / DeWitt / finite HDA structure
 -> Lorentzian finite coefficient/support controls
 -> spin-2 / TT reference sector
 -> complete six-dimensional quartic S4 TT quotient
 -> algebraic map from a future frozen six-vector to observables.
```

This chain is a **candidate mathematical/computational architecture**.  The physicalization chain in Section 10 is stronger and is not closed.

---

## 2. Binary and dimensional structure

For the frozen q=2 refinement route,

```text
N_g = (4*8^g + 10)/7

d_g = 3 + log2(1 - 35/(16*8^(g-1)+40))
```

and

```text
d_g < 3,
d_(g+1) > d_g,
d_g -> 3.
```

Finite registered diagnostics include approximately

```text
d_H          = 2.999229782
d_s(slice)   = 3.004393867
z            = 0.998281156
d_s(history) = 4.004393867
```

These are internal construction diagnostics.  They are not blind external predictions.

---

## 3. q=2 quantum-geometric carrier

The three nontrivial real Walsh characters of `Z2^2` give four exact regular-tetrahedron normals:

```text
sum_a n_a = 0
n_a.n_a   = 1
n_a.n_b   = -1/3   for a != b.
```

The four face spin-1/2 carriers contain a two-dimensional Gauss-singlet geometry sector.  The logical shape directions are represented by `X,Z`; logical `Y` is the orientation pseudoscalar.

On the selected 16-cell PL completion the same carriers glue across all 32 shared faces with the registered orientation/flux cancellation checks.

A representation-theoretic distinction is essential:

```text
four active q=2 states alone:       (2,1) + (1,2)
active states + no-link singlet:    (2,2) + (1,1)
```

under the registered endpoint structure.  The no-link state is therefore not cosmetic; it supplies the exact graph-changing completion used by the q=2 graph-link factorization.

---

## 4. Metric and Einstein reference sector

Independent finite/exact gates implement

```text
face qubits
 -> B field
 -> simplicity
 -> Urbantke metric
 -> compatible connection
 -> curvature
```

and include a non-Einstein negative control.

The logical `X/Z` shape doublet has an exact rank-two trace-free metric Jacobian.

The L1 q4 S4 metric compression gives

```text
lambda_E          = 1.1111917875584736
lambda_T2         = 1.0220278507464782
Delta_ET          = 0.08916393681199541
relative_ET_split = 0.08359564595312347
```

The ~8.36% quantity is a finite Euclidean tangent-sector split under the declared normalization.  It is not automatically a physical Lorentz-violation coefficient, particle mass ratio or measured observable.

---

## 5. Peter-Weyl constraint dynamics

The finite quantum-geometry stack contains graph/spin-changing Euclidean Hamiltonian action, volume/extrinsic-curvature controls, parity/support identities, full 32D master normalization, j=1 S4 blocking and the completed 32D higher-shell matrix.

Registered higher-shell values include

```text
lambda_min(Lambda) = 10.635759878291307
lambda_max(Lambda) = 15.059927665966466
relative distance from scalar I = 0.09440461833276048
block-Lanczos residuals ~ 1e-13.
```

These are finite **constraint-dynamics** spectral/Krylov data.

They are not particle masses and they are not a physical-frequency graviton spectrum.

In particular, for the exact Feshbach object

```text
G_c(z) = Q0^dagger (z-H_constraint)^(-1) Q0
```

`z` is a constraint spectral parameter.  It must not be renamed physical `omega` without an independently derived physical-history construction.

---

## 6. ADM / HDA / Lorentzian finite structure

The structural package includes DeWitt/ADM selection, route/diffeomorphism controls, route-normal structure function, two-node and three-node graph-changing Peter-Weyl HDA scaling, exact finite-word cutoff support, and Lorentzian coefficient/support controls.

The frozen three-node diagnostic reports approximately

```text
route exponent          = 0.9999571195
cross exponent          = 1.0024037289
pure-geometry exponent  = 2.0061524985
joint exponent          = 1.0064429344
joint defect @ 1/64     = 0.02522380790
minimum graph-change fraction = 0.4440331635.
```

This is meaningful finite off-shell evidence, but it is not a theorem uniform over arbitrary graph families, arbitrary habitats, all refinement levels and the complete Lorentzian continuum domain.

HDA consistency is also not a substitute for choosing the physical inner product or history measure.

---

## 7. Regge held-out control

The preregistered continuation

```text
Z_L = 1/8 + C/L^2 + D/L^4
```

was fitted on `L=3,4,5` and predicted

```text
Z6_pred = 0.11876923193907167
```

before comparison with

```text
Z6_obs = 0.11876075461190198
relative error ~ 0.00714%.
```

This is a genuine internal held-out numerical control.  It remains internal: it is not a held-out observation of nature.

---

## 8. TT and quartic observable dictionary

The reduced TT reference kernel is massless and has the registered positive residue / inverse-momentum equal-time covariance controls.

This reduced kernel is a **positive control**, not the final interacting theory-specific 1PI graviton kernel.

The generic parity-even S4 quartic TT quotient has exactly six independent physical structures.  The frozen six-observable extraction system has

```text
rank(100,110,111) = 5
rank(+120)        = 6
det A             = 1/699840000.
```

Thus the algebraic dictionary

```text
frozen physical six-vector
 -> TT polarization eigenvalues
 -> modified-dispersion coefficient
 -> group velocity
 -> accumulated phase / birefringence observables
```

is executable.

The crucial qualifier is `frozen physical six-vector`: the interacting microscopic values `(c1,...,c6)_IR` have not yet been derived from a theory-specific physical TT pole.

---

## 9. Finite projector and relational positive controls

Three additional results now make the legal physicalization route executable as mathematics without pretending it has already been completed for gravity.

### 9.1 Finite master-constraint theorem

For finite regulated constraints `C_A` and every positive-definite label metric `G`, define

```text
M_G = C_A^dagger G^AB C_B.
```

Then exactly

```text
M_G >= 0
ker(M_G) = intersection_A ker(C_A).
```

When zero is isolated,

```text
P_phys^(epsilon) = 1_{0}(M_G)
```

is an exact finite orthogonal projector.  Heat-kernel convergence is controlled by the first positive master gap.

What is not proved by this theorem is existence/convergence of the candidate-theory continuum/refinement rigging map.

### 9.2 Relational-history positive control

A finite C8 Page-Wootters/rigging-map model shows that

```text
global combined-constraint invariance
```

can coexist with

```text
nontrivial conditional system evolution.
```

It uses an externally declared finite clock and `R=J` as a q=2 system positive control.  Neither is claimed to be the physical gravity history generator.

### 9.3 Metric-source positive control

On that finite physical-history toy sector the legal order

```text
P_rel
 -> gauge-invariant O_rel
 -> Z[J]
 -> W[J]
 -> connected metric response
 -> tangent Gamma^(2) pseudoinverse
```

is exact.

Again, this finite `Gamma^(2)` is not the spacetime 1PI graviton kernel.

---

## 10. The physicalization chain that remains open

The no-shortcut physical route is

```text
actual regulated graph-changing constraints {C_A^(epsilon)}
 -> positive master constraint M_epsilon
 -> finite zero-sector projector
 -> controlled refinement / rigging-map or boundary-history limit
 -> frozen physical boundary/semi-classical state prescription
 -> metric-source Z[J_g]
 -> connected W[J_g]
 -> Gamma[g]
 -> physical Gamma^(2)_metric(omega,k)
 -> TT projection
 -> physical K_TT(omega,k)
 -> frozen (c1,...,c6)_IR
 -> one common physical scale
 -> preregistered blind external comparison.
```

Current machine truth is

```text
PHYSICAL_PROJECTOR_HISTORY   = open_physical
CONNECTED_INTERBLOCK_HISTORY = open_physical
PHYSICAL_TT_KERNEL           = open_physical
IR_SIX_VECTOR                = open_physical
COMMON_SCALE_CALIBRATION     = open_physical.
```

These are real physicalization gates, not cosmetic extensions.

---

## 11. Stronger extensions versus physical blockers

An arbitrary-graph HDA theorem, unbounded refinement theorem, broad microscopic universality theorem or unique blocking-measure derivation would strengthen generality.

They are different from the physical blockers in Section 10.

A finite structural candidate may be internally assembled while its theory-specific physical projector/history remains open.  Conversely, proving arbitrary-graph universality would still not automatically produce a physical `omega`-dependent 1PI kernel.

---

## 12. Experimental boundary

The repository does **not** claim

```text
experimental confirmation of quantum gravity
an observed physical graviton dispersion correction
physical cosmological constant Lambda=3
particle masses from the current gravity-only package
Standard-Model matter closure
a Born-rule derivation from the present finite positive controls.
```

External tests must use frozen outputs, one common scale rule and no post-hoc retuning.

---

## 13. What GREEN means

The canonical structural workflow `.github/workflows/core-regression.yml` answers:

```text
did the registered structural candidate reproduce its exact/finite/conditional certificates?
```

The independent `.github/workflows/physicalization-truth.yml` answers:

```text
do the solved finite projector/relational controls pass,
AND does the repository still truthfully expose the stronger physical gates as open?
```

Therefore the correct interpretation is

```text
core-regression GREEN
    = structural internal candidate reproduced

physicalization-truth GREEN
    = projector/relational reference controls reproduced
      + open physical frontier represented truthfully

both GREEN
    != experimental truth established.
```

This distinction is part of the scientific result, not merely repository bookkeeping.
