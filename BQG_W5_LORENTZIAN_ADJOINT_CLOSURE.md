# W5-seeded Lorentzian adjoint-closure runner

Date: 2026-11-10

## Purpose

The v09 adapter computes five W5 input columns but does not close its output space under the raw Lorentzian map or its adjoint. The v10 runner constructs the *reachable finite component* generated from the exact W5 seed set
`((4,)*10, (K,0,0,0,0))`, K in {0,2,4,6,8}, for one selected source node.

At each basis state it computes both:

- the raw epsilon-oriented forward map (B_v), using the action order (C_c(V)), then (C_b(K)), then (C_a(K));
- the independently ordered adjoint (B_v^dagger), using (C_a(K)), then (C_b(K)), then (C_c(V)), with the auxiliary indices transposed.

The two anti-Hermitian K-leg signs cancel. The runner does not use cyclicity of the partial auxiliary trace.

## Closure algorithm

Starting from the five W5 seed states, breadth-first search adds every output key with amplitude above the declared `--amp-tol` from either (B_v) or the direct reversed-order map. It processes each discovered key once. If `--max-states` or `--max-depth` is reached, the result is explicitly `INCOMPLETE_CLOSURE_LIMIT_REACHED`; no adjoint verdict is returned.

When closure completes, all columns are reordered into one deterministic orthonormal Gauss-key basis. The runner compares the direct reversed-order columns with the conjugate transpose of the assembled forward matrix:
[
r_{adj}=\frac{\|B_{direct\ adjoint}-B_{matrix}^{\dagger}\|_F}
{\max(\|B_{matrix}\|_F,10^{-30})}.
]
It also reports the norms of the Hermitian and anti-Hermitian parts without choosing either as the physical Lorentzian operator.

## Exact implementation details

- The v09 epsilon convention is retained: all 24 signed ordered triples.
- The direct adjoint reverses the operator order and swaps the auxiliary indices:
  `K_a[j,i] -> K_b[k,j] -> V_c[i,k]`.
- The scalar-channel restriction (J=0,1) is applied only after the middle leg, because the final covariant leg has rank 0 or 1.
- Covariant outputs are mapped back to Gauss keys only after the auxiliary trace. Rejected nonscalar outputs are counted and cause the gate to fail.
- The existing sine-ordering module installs the exact (K_{sine}=[V,H_E^{sine}]) implementation once around the full closure walk, preserving its local caches during the walk.
- Component cutoff leakage is aggregated across every processed state and both directions.

## Invocation

A bounded pilot:

```bash
python scripts/bqg_v10_w5_lorentzian_adjoint_closure.py \
  --node 0 --Jmax2 10 --max-states 40 --max-depth 2 \
  --output verification_results/W5_ADJOINT_CLOSURE_NODE0_J10_PILOT.json
```

A larger run can raise the state and depth limits. Cutoff convergence must use separate runs at `Jmax2=10,12,14`, with the same amplitude tolerance and sufficiently large closure limits. A result is cutoff-comparable only when each run reports a complete closure and the omitted-support count is zero.

## What a PASS does and does not mean

A PASS means that, on the computed W5-seeded reachable component at the declared finite cutoff, the directly assembled reversed-order map agrees with the adjoint of the assembled raw KKV matrix to the reported tolerance, and the recorded scalar-projection and component-leakage gates pass.

It does **not** establish:
- that the chosen raw KKV operator is already the physical Hermitian Lorentzian constraint;
- which symmetric ordering should be used in the full Hamiltonian;
- the (eta,kappa,hbar), regulator, or trace-normalization prefactor;
- cutoff convergence or a continuum limit;
- a physical joint kernel of (H_E+H_L).

The unit tests validate the operator-index algebra and sparse matrix assembly on finite random matrices. They do not execute the expensive Peter-Weyl closure on W5.
