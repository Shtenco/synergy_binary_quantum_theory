# BQG physical-history adapter — shortest no-new-formalism route

Status: **executable adapter contract.** This document does not redefine the physicalization architecture. It connects the actual BQG graph-changing constraint family to the already existing relational Hilbert/history, source, Gamma and observable machinery.

## 1. What is already done and must be reused

The repository already contains:

1. the finite master-constraint theorem and spectral physical projector;
2. an exact C8 constrained relational Hilbert/history sector;
3. exact gauge-invariant relational source operators;
4. an executable `P -> O_rel -> Z[J] -> W[J] -> Gamma^(2)` finite source bridge;
5. a geometric-boundary-time protocol for physical `omega`;
6. TT, six-Wilson and real-observable dictionaries;
7. Regge/TT normalization controls and a reduced connected TT propagator positive control;
8. tested-finite background `rho -> w`, scalar `mu/Sigma/slip`, and lensing/Fermat/interference reference maps.

Therefore none of these objects is to be re-invented for scalar cosmology or dark-sector work.

## 2. The only microscopic adapter that remains

On one common regulated habitat define the actual total node constraints

\[
\boxed{
H_v^{\rm BQG}=H_v^E+\lambda_L H_v^L
}
\]

with the Lorentzian ordering, epsilon convention and relative coefficient frozen before looking at the physical spectrum. The exact master is

\[
\boxed{
\mathbb M_H=\sum_v (H_v^{\rm BQG})^\dagger H_v^{\rm BQG}.
}
\]

### Exact-zero sector and the diffeomorphism generator

If the tangential/diffeomorphism action is the HDA generator derived from the same node constraints, then for an exact common zero state

\[
H_v|\psi\rangle=0\quad\forall v
\]

one has identically

\[
[H_v,H_w]|\psi\rangle=0\quad\forall v,w.
\]

Hence a commutator-derived `D` block cannot further shrink the exact common kernel. In the exact-zero case `D` is therefore an **HDA/anomaly validation**, not a second copy of the same kernel filter.

This does **not** license dropping the diffeomorphism test in a near-zero/refinement construction. If only

\[
\lambda_{low}(\mathbb M_H)\to0,
\]

then direct convergence such as

\[
\boxed{\|D\,P_{low}\|\to0}
\]

(or an equivalent uniform operator-norm theorem) is additionally required. `scripts/hda_kernel_redundancy_gate.py` contains both an exact positive control and a near-zero counterexample where the Hamiltonian master residual vanishes but a commutator residual stays finite.

The exact or controlled near-zero physical sector is then

\[
\boxed{P_{BQG}}.
\]

No dark-sector, lensing or cosmological parameter is used in constructing this projector.

## 3. Reuse the existing history/source stack unchanged

Once `P_BQG` is available, the downstream path is the already registered one:

\[
\boxed{
P_{BQG}
\to O_{rel}
\to Z_{BQG}[J]
\to W_{BQG}[J]
\to \Gamma_{BQG}.
}
\]

`scripts/bqg_physical_history_adapter_gate.py` regression-tests this interface on the existing C8 relational control by reconstructing its physical projector from the constraint `C=I-G` alone. The reconstructed projector agrees with the analytic `P_rel` to machine precision, and the resulting source and metric Hessians reproduce the existing exact q=2 source gate.

The production replacement is therefore

```text
C8 combined constraint
        ->
actual common-habitat BQG H_v = H_E,v + lambda_L H_L,v
```

as a change of microscopic input, not a change of physicalization formalism.

## 4. Three sectors of one generating functional

Dark matter, lensing and dark energy are not separate phenomenological models. They are possible sectors of the same `Gamma_BQG`.

### Scalar finite-frequency/momentum sector

After the already required scalar gauge/Dirac reduction,

\[
\Gamma^{(2)}_{scalar}(\omega,k)
\]

classifies the result as:

- analytic/contact response;
- constraint/Poisson long-range response;
- or an additional genuine propagating scalar pole.

Only the last case can be discussed as an extra physical scalar degree of freedom. A modified Poisson response can mimic dark-matter phenomenology without a new particle, but must also pass the universal conserved-source and lensing tests.

The existing `physical_cosmology_background_scalar_gate.py` is reused only as the downstream conservation and `mu/Sigma/slip` dictionary; it does not supply the BQG scalar kernel.

### TT sector

\[
\Gamma^{(2)}_{TT}(\omega,k)
\]

feeds the already closed six-Wilson and real-observable map. No new TT observable dictionary is needed.

### Homogeneous background sector

\[
\Gamma_{FLRW}[a,N]
\]

must come from absolute physical history amplitude ratios with one fixed normalization convention. Background-by-background normalization is forbidden because it can erase a genuine volume/vacuum term.

From the same functional one derives

\[
\rho_{hist}(a),\quad p_{hist}(a),\quad w_{hist}(a).
\]

The existing background gate is intentionally one-way: `rho_hist(a) -> p_hist(a), w_hist(a)`. A desired `w(a)` is never used to choose `rho_hist`.

### Lensing / interference consistency

The existing `binary_history_interference_lensing_gate.py` freezes the rule that the **same** Weyl/Fermat potential must control deflection, stationary images, time delay and coherent phase. Its split-potential example is a negative control. It is not a theory-specific BQG lensing result.

## 5. Production bundle contract

Every actual BQG constraint-history bundle must freeze before diagonalization:

- regulator/refinement level;
- common habitat/basis identity;
- `H_E` ordering and normalization;
- `H_L` ordering, epsilon convention and relative coefficient convention;
- HDA/diffeomorphism validation convention;
- boundary embedding and refinement map;
- source-operator definitions;
- hashes of all input matrices/operators.

The adapter must fail closed if a required microscopic component is missing or acts on a different habitat.

## 6. Forbidden shortcuts

The following do not count as physicalization:

```text
P_boundary M P_boundary -> call its kernel P_phys
constraint resolvent z -> rename z as omega
K1 Gram eigenvalue -> call it physical mass/gap
j=1 volume eigenvalue -> call it dark energy
observer smoothing b^-2 -> call it 1/k^2 propagator
missing H_L -> fit a phenomenological correction
near-zero Hamiltonian mode -> assume D/HDA residual is also zero
normalize Z[a;0]=1 separately for every background -> infer Lambda=0
fit a separate optical/lensing potential after seeing dynamics
```

## 7. Current actual BQG input status

The regulator-safe Euclidean calculations show that the bare 32D q=2 boundary carrier has no common five-node Euclidean zero direction and that its first outgoing layer has full rank 160. These are input diagnostics, not physical observables.

The first 8-shard full 24-term Lorentzian epsilon run hit its 90-minute runner limit. The scientific protocol is unchanged and has been rescheduled as 24 one-term shards with the same 24 signed ordered terms and the same `Jmax=7/2`.

The Euclidean `K1` Ritz second-action calculation is running independently. Its result is a representation-dressing diagnostic, not yet `P_BQG` because the common-habitat total Lorentzian node constraints remain required.

## 8. Shortest legal programme

\[
\boxed{
\{H_v^E+\lambda_L H_v^L\}
\xrightarrow[\text{HDA validation}]{\text{exact/controlled kernel}}
P_{BQG}
\to \text{existing relational/history stack}
\to \text{existing }Z/W/\Gamma
\to
\begin{cases}
\Gamma^{(2)}_{scalar}(\omega,k),\\
\Gamma^{(2)}_{TT}(\omega,k),\\
\Gamma_{FLRW}[a,N]
\end{cases}
\to \text{existing cosmology/lensing/TT/real-observable dictionaries}.
}
\]

This is the production programme. No dark-matter or dark-energy parameter is allowed to enter before the final sector diagnostics are computed.
