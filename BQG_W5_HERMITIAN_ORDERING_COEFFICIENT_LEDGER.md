# W5 Lorentzian Hermitian-ordering and coefficient ledger

**Date:** 2026-11-11  
**Scope:** exact algebraic conversion of the current raw W5 KKV map; comparison with a published real-Ashtekar–Barbero regularization; no physical prefactor is declared fitted or finalized.

## 1. Frozen operator actually used by V10

The current W5 V10 runner uses
\[
H_E^{sine}=\frac{T-T^\dagger}{2i},\qquad
K_{comm}=[V,H_E^{sine}]=VH_E^{sine}-H_E^{sine}V.
\]
On a common finite domain, with \(V^\dagger=V\) and \((H_E^{sine})^\dagger=H_E^{sine}\),
\[
K_{comm}^\dagger=-K_{comm}.
\]
The raw node-summed map is
\[
B_{comm}=\sum_v\sum_{(a,b,c)\in\mathcal O_v}
\epsilon(a,b,c)\operatorname{Tr}_{aux}
[C_a(K_{comm})C_b(K_{comm})C_c(V)].
\]
The direct reversed-order map computed by V10 is denoted \(B_{comm}^{\dagger,dir}\). Its matrix must agree with the conjugate transpose of the fully assembled \(B_{comm}\) on one common adjoint-closed basis before any ordering completion is accepted.

## 2. Exact Hermitian K conversion and the sign

If the commutator is the quantum representation of the corresponding Poisson bracket in the adopted convention, the Hermitian K-leg candidate is
\[
K_H=\frac{1}{i\hbar}K_{comm}.
\]
Then
\[
K_{H,a}K_{H,b}V_c
=\frac{1}{(i\hbar)^2}K_{a,comm}K_{b,comm}V_c
=-\frac{1}{\hbar^2}K_{a,comm}K_{b,comm}V_c.
\]
Thus the exact algebraic conversion of the two K legs is
\[
B_H=-\hbar^{-2}B_{comm},
\]
provided no additional operator factors or ordering changes are introduced. The minus sign is mandatory. It remains at \(\hbar=1\).

This identity alone does **not** prove that the repository's \(K_{comm}\) has the complete canonical normalization of extrinsic curvature. Its sign relative to the classical Poisson-bracket convention, and its missing \(\kappa,\beta\), regulator and volume factors, must be checked against one frozen canonical convention.

## 3. Symmetric completion of the chosen raw ordering

Once the full adjoint-closed matrix is available, the minimal symmetric completion of this specific raw ordering is
\[
B_{comm,sym}=\frac12(B_{comm}+B_{comm}^\dagger),
\qquad
B_{comm,anti}=\frac12(B_{comm}-B_{comm}^\dagger).
\]
The first is Hermitian; the second is anti-Hermitian. Both must be retained in the audit. The corresponding two-Hermitian-K candidate is
\[
B_{H,sym}=-\frac{1}{2\hbar^2}
(B_{comm}+B_{comm}^\dagger)
=-\hbar^{-2}B_{comm,sym}.
\]

This is a mathematically unambiguous symmetric completion of the selected raw operator. It is **not yet proven to be the unique physical Thiemann ordering** for the BQG discretization. A different placement of holonomies, volume factors, inverse-volume factors, or node-lapse weights can define a different symmetric operator with different quantum corrections.

## 4. Primary literature comparator and why coefficients cannot be copied

Yang and Ma, *New Hamiltonian constraint operator for loop quantum gravity*, Phys. Lett. B 751 (2015) 343–347, DOI 10.1016/j.physletb.2015.10.062, write the real Ashtekar–Barbero Hamiltonian as \(H=H_E-T\) with a Lorentzian term carrying \(1+\beta^2\). Their regulated Lorentzian term (their Eq. 9) contains convention-specific factors and inverse-volume insertions. Their operator construction (Eqs. 12–13) defines the extrinsic-curvature operator through a separately normalized commutator involving \(H_E\), \(V\), \(\beta\), and Planck-scale/representation factors. They then define the symmetric regulated Hamiltonian by
\[
\hat H_{sym}=\frac12(\hat H+\hat H^\dagger)
\]
(their Eq. 20; the isotropic model is symmetrized again in Eq. 26).

This gives a legitimate **methodological comparator**: the full regulated Hamiltonian is symmetrized after its own Lorentzian operator and normalization have been constructed. It does not authorize copying their numerical coefficient into the BQG KKV matrix. Their operator contains a particular regularization, inverse-volume factors, holonomy placement, representation trace and definition of \bar K, none of which has yet been shown identical to the present finite KKV map.

Primary sources:
- Journal/publisher record: https://www.sciencedirect.com/science/article/pii/S0370269315008217
- Open full text: https://scoap3-prod-backend.s3.cern.ch/media/files/12442/10.1016/j.physletb.2015.10.062.pdf
- arXiv: https://arxiv.org/abs/1507.00986

## 5. Classical relative coefficient already fixed in this repository

The repository's classical real-\(\beta\) kinetic-sector identity is
\[
H_E^{kin}=-\beta^2Q_{DW},\qquad
H_L^{corr}=(1+\beta^2)Q_{DW},
\qquad
H_E^{kin}+H_L^{corr}=Q_{DW}.
\]
This fixes the relative classical coefficient in the derivative-free homogeneous kinetic sector. It does not determine the absolute coefficient of the finite W5 KKV matrix elements.

The only safe ledger form at present is therefore
\[
\hat H_{L,W5}
=\mathcal N_{W5}^{reg}
(\beta,\kappa,\hbar,\ell_P,\ell,\text{trace},\text{orientation},N)
(1+\beta^2)B_{H,sym},
\]
where \(\mathcal N_{W5}^{reg}\) is **unknown and must be derived**, not fitted. Depending on the frozen operator definition, some factors may be absorbed into \(B_{H,sym}\); the ledger must prevent double-counting.

## 6. Required convention ledger before a physical joint kernel

The production artifact must record, separately:

| Item | Current status |
|---|---|
| \(H_E^{sine}=(T-T^\dagger)/(2i)\) | Implemented structural ordering |
| \(K_{comm}=[V,H_E^{sine}]\) | Implemented; anti-Hermitian |
| Direct reversed-order \(B_{comm}^\dagger\) | Implemented in V10; production closure pending |
| \(B_{comm,sym}\) and \(B_{comm,anti}\) | Algebraically defined; full W5 numerical matrices pending |
| \(K_H=K_{comm}/(i\hbar)\) | Exact algebraic conversion conditional on canonical commutator convention |
| \(B_{H,sym}=-B_{comm,sym}/\hbar^2\) | Exact two-leg conversion under the same condition |
| Classical relative \(1+\beta^2\) | Fixed in the repository's kinetic-sector identity |
| \(\kappa\), Planck-scale, representation-trace and inverse-volume factors | Not mapped to the current KKV normalization |
| Edge/node orientation and lapse weights | Must be frozen in one production convention |
| Regulator/triangulation-cell conversion | Not established |
| Unique physical BQG ordering | Open |

No coefficient may be adjusted to force a joint kernel, maximize a small singular value, or reproduce a target spectrum.

## 7. Acceptance sequence

1. Complete V10 closure and direct-adjoint comparison at \(J_{max,2}=10,12,14\), with identical seed basis, amplitude tolerance, node set and closure limits.
2. Require complete closure at each cutoff; otherwise compare only explicitly labelled partial diagnostics.
3. Report the full matrices \(B_{comm}\), \(B_{comm}^\dagger\), symmetric and anti-Hermitian parts, singular spectra, rank, condition numbers and leakage.
4. Independently verify the operator-domain and adjoint identities of every K/V leg, not just the final matrix.
5. Derive \(\mathcal N_{W5}^{reg}\) by translating one selected primary classical regularization into the exact code conventions. Do not transplant a literature prefactor without mapping every factor.
6. Only then form the common-carrier physical \(H_E+H_L\) and its joint kernel.

## 7.1. V10 implementation update (11 November 2026)

The V10 runner has now been extended to serialize, on a complete adjoint-closed basis, three separate sparse column families:
- (B_{comm,sym}=(B_{comm}+B_{comm}^\dagger)/2);
- (B_{comm,anti}=(B_{comm}-B_{comm}^\dagger)/2);
- (B_{H,sym}|_{\hbar=1}=-B_{comm,sym}).

The runner explicitly labels the last object as a dimensionless (\hbar=1) algebraic candidate. It does not mark the physical ordering as selected and does not insert the missing external regularization coefficient.

A new unit test constructs a complex matrix and checks (S=S^\dagger), (A=-A^\dagger), (S+A=B), and (B_{H,sym}|_{\hbar=1}=-S). This is a matrix-algebra regression test, not a W5 production result. The corresponding Actions run must pass before this new code path is treated as tested.

## 8. Current verdict

**Established algebraically:** the minus \(\hbar^{-2}\) conversion for two anti-Hermitian commutator legs, and the unique minimal symmetric completion of the chosen raw map once its full direct adjoint is known.

**Not established:** that this is the unique physical BQG Lorentzian ordering; the absolute \(\kappa,\beta,\hbar,\ell_P\), trace and regulator normalization; the physical W5 Lorentzian matrix; the physical \(H_E+H_L\) kernel.

A successful finite adjoint-closure run is necessary structural evidence. It is not a physical Hamiltonian certificate.
