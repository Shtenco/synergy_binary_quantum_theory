# q=2 Walsh local QFI control — result (2026-11-18)

## Exact input

The existing `micro_walsh_qgeom_gate.py` constructs four pure face-qubit
states from the four route labels (00,01,10,11). Their Bloch vectors are the
regular tetrahedral frame
[
n_a\cdot n_b=-\frac13\quad(a\ne b),\qquad
\sum_a n_a=0,\qquad
\sum_a n_a n_a^T=\frac43 I_3.
]

For a pure qubit with state (|n_a\rangle) and SU(2) frame generators
(J_i=\sigma_i/2), the unitary-parameter QFI matrix is
[
(F_Q^{(a)})_{ij}
=4\operatorname{Re}\bigl(
\langle J_iJ_j\rangle-\langle J_i\rangle\langle J_j\rangle
\bigr)
=\delta_{ij}-n_{a,i}n_{a,j}.
]

The exact four-face sum is therefore
[
\sum_{a=1}^4F_Q^{(a)}
=4I_3-\sum_a n_an_a^T
=\boxed{\frac83 I_3}.
]
The average per face is (\frac23I_3). This is a genuine local
information-geometric isotropy check derived from the existing Walsh states.

## Gauss projection changes the answer

The canonical four-face product state has Gauss-singlet projection weight
[
\langle\psi|P_{\rm singlet}|\psi\rangle=\frac29.
]

After projection and normalization, the state is annihilated by all three
global generators
[
J_i^{\rm tot}=\frac12\sum_{f=1}^4\sigma_i^{(f)}.
]
The numerical residual is (1.12\times10^{-16}), and the global gauge-orbit
QFI matrix is zero up to floating-point roundoff ((<4\times10^{-17}) per
diagonal entry).

**Consequence:** a global SU(2) gauge rotation cannot be used as a positive
QFI edge conductance after Gauss projection. A nonzero network metric must
come from a relational/inter-block generator, not from the gauge orbit itself.

## Reduced relative-twist control

For the unprojected product face states, the reduced-model probe
[
K_{ij}=\frac12(Z_i-Z_j)
]
gives, for all six pairs,
[
F_Q(|\psi\rangle,K_{ij})=\frac43,
\qquad
g_{ij}^{\rm probe}=\frac14F_Q=\boxed{\frac13}.
]

All six values agree to numerical precision. This is a useful positive
control, but it is not a physical edge conductance: the state is a product
state, the generator is a reduced-model relative twist, and the pairs are
local face-qubit pairs rather than derived inter-block physical links.

Moreover, a common scalar (g=1/3) only rescales the graph diffusion
parameter. It cannot select spectral dimension. Any (d_s) obtained by
putting this weight on a chosen graph would still be a property of that
graph's topology and measure.

## CI

GitHub Actions run `44271370461` passed:
- 3 tests;
- JSON result export;
- artifact upload `q2-walsh-local-qfi-control`.

The result is reproducible from:

```bash
python scripts/q2_walsh_local_qfi_control.py \
  --output verification_results/Q2_WALSH_LOCAL_QFI_CONTROL.json
```

## Claim boundary and next requirement

Closed: local Walsh-frame QFI isotropy; zero global gauge QFI after the
Gauss-singlet projection; exact reduced product-state relative-twist control.

Open: physical two-block state (ρ_{e}), relational generator (K_e),
gauge-covariant gluing across a shared tetrahedral face, and a derivation of
edge conductance from the physical constraint/history amplitude. Do not use
(1/3) as a physical QFI weight on the PL dual graph.

The current PL exporter now also exports the tetrahedral dual graphs, with
route-coloured shared-face edges:
- generation 0: 16 dual vertices, 32 edges;
- generation 1: 384 dual vertices, 768 edges;
- generation 2: 9,216 dual vertices, 18,432 edges.

Every dual graph is connected and 4-regular; each shared face has exactly two
incident tetrahedra with matching route colour. This is the natural PL
candidate adjacency for a future relational QFI calculation, but it remains
distinct from the scalar causal rewrite DAG.
