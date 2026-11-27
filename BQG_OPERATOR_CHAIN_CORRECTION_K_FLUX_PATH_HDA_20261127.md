# BQG operator-chain correction: existing K block, flux sharp map, and path-normal HDA — 2026-11-27

Status: **corrective audit and integration plan**. This document supersedes the overly broad sentence in `BQG_NORMAL_EVOLUTION_EXTRINSIC_CURVATURE_OPERATOR_AUDIT_20261127.md) that said no microscopic extrinsic-curvature operator had been found. A finite commutator-derived K-leg block does exist. It is a prerequisite block, not yet the spatial tensor \(K_{ab}\) or a physical Lorentzian normal-evolution generator.

## 1. Correction to the previous audit

The inspected `main` branch contains the following concrete constructions:

| Component | Exact object | Existing evidence | Boundary |
|---|---|---|---|
| K-leg | \(K_{v,comm}=[V_v,H^{sine}_{E,v}]\) | `scripts/peter_weyl_lorentzian_K_block_gate.py`; W5 adapter uses the same sine-K engine | Local finite commutator; not yet the full tensor \(K_{ab}\) |
| Hermitian K conversion | \(K_{v,H}=(i\hbar)^{-1}K_{v,comm}\) | `BQG_W5_K_COMMUTATOR_NORMALIZATION.md` | Exact conversion, but absolute canonical prefactors remain open |
| Lorentzian KKV map | \(\sum\epsilon\,\mathrm{Tr}_{aux}(K_aK_bV_c)\) | `scripts/bqg_v10_w5_lorentzian_adjoint_closure.py` | Raw structural map; ordering, normalization, scalar covariance, and physical joint Gram not closed |
| Lapse cochain | \(\omega_{vw}=N_vM_w-N_wM_v\) | `scripts/dual_k5_lapse_cochain_gate.py` | Exact combinatorics; no metric raising by itself |
| Flux/RT0 reconstruction | dual one-cochain \(\to\) tetrahedral vector field | `scripts/dual_cell_sharp_rt0_gate.py` | Exact local kinematics for circumcentric duals; not yet the full state-dependent quantum HDA map |
| Path diffeomorphism | \(\hat D_{path}[\beta]\) with two transverse route directions | `scripts/path_vector_diffeo_gate.py`, `scripts/path_rerouting_diffeo_gate.py` | Kinematic route representation; not yet coupled to the same quantum geometry/constraint domain |
| Path-normal generator | \(H_{path}[N]=\tfrac12\{N,\sqrt{-\Delta_{path}}\}\) | `scripts/path_normal_hda_gate.py` | Correct principal-symbol HDA route-sector construction; not the full Peter–Weyl gravity Hamiltonian |

The previous audit should therefore be read as: **no fully assembled, physical, Lorentzian normal-evolution generator and no derived spatial \(K_{ab}\) have yet been certified**. It should not be read as saying that no K-related operator exists.

## 2. Exact K-commutator bookkeeping

The current sine-Euclidean operator is
\[
H_E^{sine}=\frac{T-T^\dagger}{2i},
\qquad
K_{comm}=[V,H_E^{sine}].
\]
Both entries in the commutator are Hermitian, hence
\[
K_{comm}^\dagger=-K_{comm}.
\]
The Hermitian operator corresponding to this commutator is
\[
\boxed{K_H=\frac{1}{i\hbar}K_{comm}=-\frac{i}{\hbar}K_{comm}.}
\]
For two K legs,
\[
\boxed{K_{H,a}K_{H,b}V_c=-\hbar^{-2}K_{comm,a}K_{comm,b}V_c.}
\]
The minus sign is mandatory. It cannot be absorbed into an undocumented Lorentzian coefficient or fitted from the Gram spectrum.

This establishes an exact algebraic conversion of the implemented commutator. It does **not** yet establish that the finite operator has the complete canonical normalization of \(K_a^i\), nor that a spatial tensor \(K_{ab}\) has been reconstructed. The latter needs a triad/flux-to-tensor map, a common frame/transport convention, and the physical-state/continuum limit.

## 3. W5 is still blocked before a physical KKV result

The W5 runner uses the correct five spin-2 seed states
\[
((4)^{10},(K,0,0,0,0)),\qquad K\in\{0,2,4,6,8\},
\]
and constructs forward KKV columns plus independently reversed-order adjoint columns. Its closure walk is designed to retain the actual amplitudes and reject non-scalar output.

However, the W5 scalar-covariance pilot recorded a small nonzero rejected-sector weight:
- expected \(J=0,M=0\) weight: \(6.641081550043974\);
- rejected \(J=1,M=0\) weight: \(5.390512613173637\times10^{-20}\);
- relative amplitude scale: approximately \(9.0094\times10^{-11}\).

The existing gate correctly refuses to silently project this away. The cause has not been classified as floating-point cancellation, recoupling/orientation mismatch, or a real finite-cutoff covariance defect. Until the independent oriented-spec, T/T-dagger, cutoff, and precision decomposition is complete, the W5 Lorentzian KKV output is not certified.

The W5 joint-Gram runner also correctly requires raw, row-aligned complex \(A\) and \(B\) columns. The cross-Gram \(A^\dagger B\) cannot be reconstructed from separate Gram matrices. No physical W5 joint kernel is claimed.

## 4. The path-normal HDA construction is a second, independent route

The path-normal script defines
\[
\Omega=\sqrt{-\Delta_{path}},\qquad
H_{path}[N]=\frac12\{N,\Omega\}.
\]
Its principal symbol is
\[
h_N=N\sqrt{q^{ab}p_ap_b}.
\]
The Poisson bracket of these symbols gives the standard metric structure function (up to one fixed global vector-constraint orientation convention):
\[
\{h_N,h_M\}
=q^{ab}(M\partial_bN-N\partial_bM)p_a.
\]

The numerical gate compares the finite spectral commutator against the corresponding vector generator on WKB carriers and checks the expected high-frequency defect decay. This is meaningful: it provides a constructive normal-deformation representation in the route sector without fitting the magnitude of the structure function.

But the script's own scope note is decisive: the full Peter–Weyl Lorentzian geometry operator must still couple to this same path domain. Therefore this is not yet a BQG physical Hamiltonian or a derived spacetime light cone.

## 5. Flux-driven shift: what is already closed and what remains

The kinematic chain is
\[
(N,M)\to\omega_{vw}=N_vM_w-N_wM_v
\to\beta_v\to\beta_{\perp,e}\to D_{path}[\beta].
\]
For a circumcentric tetrahedral dual, the local flux reconstruction can be expressed through the face fluxes and RT0 basis. The RT0 gate tests:
- exact reconstruction of constant vector fields from dual-edge integrals;
- exact reproduction of arbitrary face fluxes;
- closure of the constant-vector face fluxes;
- orthogonality of adjacent circumcenters to the shared face.

The `FLUX_DRIVEN_PATH_HDA_TARGET.md` further proposes a cell-centred matrix
\[
G_v=\sum_{w\sim v}\frac{d_{vw}}{A_{vw}}E_{vw}E_{vw}^{T},
\qquad
\beta_v=G_v^{-1}\sum_{w\sim v}\omega_{vw}E_{vw}.
\]
This formula is a candidate state-dependent sharp map. Its acceptance must be checked against the RT0 reconstruction on the same oriented tetrahedral data, including dimensions, flux closure, neighboring-cell averaging, degeneracy conditioning, and gauge/frame covariance. The RT0 constant-field test alone does not certify this entire quantum-state-dependent formula.

## 6. The missing integration theorem

The key integration target is not simply to place the three existing scripts beside one another. They currently use different mathematical domains and levels of description. A valid closure must construct a commuting diagram in one common regulated setting:

\[
\begin{array}{ccc}
(N,M,\Psi_{PW}) &\xrightarrow{\ [\hat H[N],\hat H[M]]\ }& \mathcal H_{PW}^{reachable}\\
\downarrow\omega\;\downarrow\mathrm{flux/metric}&&\downarrow\mathrm{coarse\ graining}\\
\beta=\sharp_{E,q}\omega&\xrightarrow{\ \hat D_{path}[\beta]\ }&\mathcal H_{path}\otimes\mathcal H_{PW}
\end{array}
\]
with the same physical state, regulator, lapse normalization, and graph/route embedding on both paths.

Separately, the K-commutator must enter the same Lorentzian constraint that generates normal evolution. Its classical limit must reproduce
\[
K_{ab}=\frac{1}{2N}(\dot q_{ab}-\mathcal L_{\vec N}q_{ab})
\]
for the fixed \((-+++),\ n^\mu n_\mu=-1,\ K_{ab}=-q_a{}^\mu q_b{}^\nu\nabla_\mu n_\nu\) convention. The finite \([V,H_E]\) block is an ingredient for this test, not its conclusion.

## 7. Updated acceptance sequence

1. **K-leg validation:** resolve W5 scalar covariance with per-oriented-spec and T/T-dagger decomposition, then precision/cutoff sweeps. No scalar projection to force PASS.
2. **K normalization:** preserve the exact \(-\hbar^{-2}\) conversion for two Hermitian K legs and separately freeze \(\kappa,\beta\), regulator, edge/node, and trace factors from one classical convention.
3. **K tensor reconstruction:** define the flux/triad-to-\(K_{ab}\) map and prove frame/gauge covariance; do not identify the scalar node commutator with the tensor.
4. **Sharp-map validation:** test the proposed \(G_v^{-1}\sum\omega E\) formula against RT0 on identical geometry data and on neighboring-cell interfaces; report condition numbers and degenerate cells.
5. **Path/HDA coupling:** couple \(D_{path}\) and \(H_{path}[N]\) to the same Peter–Weyl geometry states and compare the full reachable-space commutator with the geometry-derived right-hand side.
6. **Lorentzian completion:** construct the Hermitian Lorentzian correction with fixed classical coefficient and ordering, then calculate raw common-basis amplitudes and the joint W5 Gram.
7. **Relational evolution:** derive the clock and relational generator from the same physical construction; do not identify the master-constraint parameter or path-register coordinate with physical time without proof.
8. **Metric and cone:** derive \(q_{ab},K_{ab},N,N^a\), assemble the Lorentzian metric, and only then feed its null generators into the existing optical Jacobi/Riccati gate.

## 8. Current scientific conclusion

The repository is further along than the previous audit stated: it contains a genuine commutator-derived K-leg, a tested local flux-to-vector reconstruction, a path-diffeomorphism representation, and a square-root path-normal HDA principal-symbol construction.

What remains open is their **joint physical identification**. In particular, the W5 scalar-covariance blocker, canonical normalization of the K operator, state-dependent sharp-map validation, coupling of path and Peter–Weyl domains, Lorentzian Hermitian completion, and relational-time construction remain mandatory.

No new null cone should be fitted from these separate components. A physical BQG cone follows only after the common-domain HDA and normal-evolution chain closes.
