# W5-seeded Lorentzian adjoint-closure runner

Date: 2026-11-10

## Purpose

The v09 adapter computes five W5 input columns but does not close its output space under the raw Lorentzian map or its adjoint. The v10 runner constructs the reachable finite component generated from the exact W5 seed set:

`((4,)*10, (K,0,0,0,0))`, with `K in {0,2,4,6,8}`.

By default it includes all five K5 source nodes and coherently sums their contributions. A node-local diagnostic can be run by selecting a subset of nodes.

## Forward map and direct adjoint

At each basis state and selected node, it computes:

- The raw epsilon-oriented forward map `B_v`, with right-to-left action `C_c(V) -> C_b(K) -> C_a(K)`.
- The independently ordered adjoint `B_v†`, with right-to-left action `C_a(K) -> C_b(K) -> C_c(V)`, and transposed auxiliary indices.

In explicit indices, the forward contraction is `Σ A_ij B_jk C_ki`; the adjoint is `Σ C_ik B_kj A_ji`. The two anti-Hermitian K-leg signs cancel. The implementation does not use cyclicity of the partial auxiliary trace.

For the default unit-lapse case it forms `B = Σ_v B_v` and `B† = Σ_v B_v†`. If a nonuniform lapse is intended, its node weights must be supplied and included in the adjoint construction; the current runner does not infer them.

## Closure algorithm

Starting from the five W5 seed states, breadth-first search adds every output key with amplitude above the declared `--amp-tol` from either the coherent forward map or its direct reversed-order adjoint. Each discovered key is processed once.

If `--max-states` or `--max-depth` is reached, the result is explicitly `INCOMPLETE_CLOSURE_LIMIT_REACHED`; no adjoint verdict is returned.

When closure completes, all columns are reordered into one deterministic Gauss-key basis. The runner compares the direct reversed-order columns with the conjugate transpose of the assembled forward matrix:

`r_adj = ||B_direct_adjoint - B_matrix†||_F / max(||B_matrix||_F, 1e-30)`.

It also reports the Frobenius norms of the Hermitian and anti-Hermitian parts without choosing either as the physical Lorentzian operator.

## Exact implementation details

- The v09 epsilon convention is retained: all 24 signed ordered triples at each source node.
- The scalar-channel restriction J=0,1 is applied only after the middle leg, because the final covariant leg has rank 0 or 1.
- Covariant outputs are mapped back to Gauss keys only after the auxiliary trace. Rejected nonscalar outputs are counted and cause the gate to fail.
- The Gauss-to-covariant embedding is a one-to-one relabeling of the source-node singlet into the J=0, M=0, K12=K34 covariant key. Its norm is preserved by construction; the scalar output map is its inverse on that sector.
- The existing sine-ordering module installs the exact K_sine=[V,H_E^sine] implementation once around the full closure walk, preserving its local caches during the walk.
- Component cutoff leakage is aggregated across every processed state, direction, and selected node.
- Amplitude pruning is deferred until the complete node-summed forward and adjoint columns have been coherently assembled.

## Invocation

A bounded all-node pilot:

```bash
python scripts/bqg_v10_w5_lorentzian_adjoint_closure.py \
  --nodes 0 1 2 3 4 --Jmax2 10 --max-states 40 --max-depth 2 \
  --output verification_results/W5_ADJOINT_CLOSURE_ALLNODES_J10_PILOT.json
```

A node-local diagnostic can pass a subset such as `--nodes 0`. A larger run can raise the state and depth limits. Cutoff convergence requires separate runs at Jmax2=10,12,14 with the same amplitude tolerance and sufficiently large closure limits. A result is cutoff-comparable only when each run reports complete closure and zero omitted support.

## What a PASS does and does not mean

A PASS means that, on the computed W5-seeded reachable component at the declared finite cutoff, the directly assembled reversed-order map agrees with the adjoint of the assembled raw KKV matrix to the reported tolerance, and the recorded scalar-projection and component-leakage gates pass.

It does **not** establish:
- that the raw KKV operator is already the physical Hermitian Lorentzian constraint;
- which symmetric ordering should be used in the full Hamiltonian;
- the beta, kappa, hbar, regulator, or trace-normalization prefactor;
- cutoff convergence or a continuum limit;
- a physical joint kernel of H_E+H_L.

The unit tests validate the operator-index algebra and sparse matrix assembly on finite random matrices. They do not execute the expensive Peter-Weyl closure on W5.
