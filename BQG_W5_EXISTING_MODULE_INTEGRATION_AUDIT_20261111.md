# W5 integration audit: reuse of existing Python modules

**Audit date:** 2026-11-11  
**Branch:** \`research/lemaitre-bqg-causal-cone-20261110\`  
**Audited HEAD before this note:** \`89eadb9873023448186508c0279e7583963687c7\`

## Executive conclusion

The W5 programme is not a blank-sheet implementation. The repository already contains three complementary computational layers:

1. V09 produces per-input, per-source-node Euclidean sine amplitudes and raw epsilon KKV Lorentzian amplitudes on W5.
2. V10 coherently sums the selected source-node raw KKV maps, closes the finite domain under both the forward map and its independently assembled adjoint, and exposes the symmetric/anti-Hermitian decomposition.
3. The joint-Gram gate combines raw Euclidean and Lorentzian maps only after they are represented in the same ordered output carrier.

The missing work is an integration and physical-normalization problem, not a reason to replace these modules with a second implementation. This audit does not claim the production amplitudes or physical kernel have been computed.

## Repository inventory actually checked

The recursive Git tree at the audited branch contained 727 tracked blobs:
- 211 \`.py\` files;
- 212 \`.md\` files;
- 77 \`.pkl\` files;
- 75 \`.json\` files;
- 50 \`.svg\` files;
- 43 \`.log\` files;
- 40 workflow \`.yml\` files;
- 8 zip files.

This is the inventory of this branch snapshot, not a claim about untracked/generated files or other repositories. Cross-checks of the repository's `main` tree found 210 tracked Python files; the historical `research/bqg-v26-multiprecision-charge-20261109` tree found 278. None of these checked Git trees exposes 1500+ tracked Python modules. The broader 1500+ inventory may refer to another workspace, generated modules, or repositories not represented by these refs; it must be located before it can be audited rather than assumed to be present here.

## Task A — W5 V10 closure and cutoff study

### Existing production path

\`scripts/bqg_v10_w5_lorentzian_adjoint_closure.py\` reuses:
- \`scripts/k5_peter_weyl_safe_hda_column.py\`;
- \`scripts/peter_weyl_covariant_volume_leg_gate.py\`;
- \`scripts/peter_weyl_lorentzian_sine_ordered_triple_gate.py\`;
- \`scripts/peter_weyl_lorentzian_ordered_triple_gate.py\`.

The seed is the five-state W5 carrier
\[
((4)^{10},(K,0,0,0,0)),\qquad K\in\{0,2,4,6,8\}.
\]
The runner coherently sums source nodes 0–4 with unit lapse, closes under forward and direct reversed-order maps, and compares the direct adjoint with the matrix adjoint on one deterministic finite basis.

### Current run state

The full-node workflow has matrix cutoffs \(J_{\max,2}=10,12,14\), with \`max-states=500\`, \`max-depth=8\`, and \`amp-tol=10^{-12}\`. The three matrix jobs in run 42897439770 were observed queued; no closure artifact was available at the time of this audit. The adapter test run 42897439803 was also queued. The additional joint-Gram regression run was triggered by the subsequent gate update; its status must be read from Actions rather than inferred from the push.

### Acceptance logic

A cutoff point is eligible for comparison only if:
- \`status == PASS_ADJOINT_CONSISTENCY\`;
- closure is complete, with zero omitted support and no depth/state cap reached;
- all five W5 seeds and all five source nodes are present;
- scalar-projection rejection is zero;
- all hard component leakages are below the declared threshold;
- the direct-adjoint residual is below tolerance;
- the basis identity, amplitude tolerance, node weights, and operator convention match the other cutoffs.

Even three eligible points do not by themselves prove convergence. Compare the same observables across 10/12/14: closure size and depth, basis nesting/overlap, forward and adjoint Frobenius norms, per-column norms and supports, component leakage, adjoint residual, symmetric/anti-Hermitian norms, numerical rank, and singular spectrum. A plateau must be stable against the amplitude-pruning threshold and increased closure caps. Incomplete runs are diagnostic only.

## Task B — physical Lorentzian amplitudes, ordering, and normalization

### Existing amplitude sources

\`scripts/bqg_v09_w5_lorentzian_raw_column.py\` returns, for one W5 input K and one source node:
- \`A_Euclidean_sine_Gauss\`;
- \`B_Lorentzian_raw_KKV_Gauss\`;
- the covariant raw output and leakage diagnostics.

The V09 result explicitly labels B as raw KKV, without the Hermitian completion and without the beta/kappa/hbar prefactor. It is therefore a source of structural amplitudes, not the physical \(H_L\).

\`scripts/bqg_v10_w5_lorentzian_adjoint_closure.py\` supplies the all-node coherent raw KKV map, its independently assembled adjoint, and
\[
B_{\rm comm,sym}=\tfrac12(B_{\rm comm}+B_{\rm comm}^\dagger),\qquad
B_{\rm comm,anti}=\tfrac12(B_{\rm comm}-B_{\rm comm}^\dagger).
\]
The current node sum is explicitly unweighted (unit lapse).

### Required integration, without rewriting the engines

The W5 common-carrier builder must:
1. run the existing V09 Euclidean sine map for each of the five W5 inputs and each of the five source nodes, retaining raw complex Gauss-key amplitudes;
2. coherently sum node contributions for each input using one declared lapse convention;
3. take the W5 seed columns of the completed V10 matrix by exact seed-key identity, not by assuming a column order;
4. use V10's direct-adjoint symmetric part for the raw Lorentzian Hermitian candidate;
5. construct one deterministic union output basis for the Euclidean and Lorentzian maps, preserving exact Gauss keys, phases, and row order;
6. record hashes for the input seeds, output basis, source modules, cutoff, amplitude tolerance, node weights, and both raw maps;
7. feed both maps into \`scripts/bqg_lorentzian_w5_joint_gram_gate.py\` only after the canonical coefficient ledger is frozen.

The raw V09 per-node B map is also an independent cross-check of the V10 forward columns. It should be compared node-by-node and after coherent summation before the Hermitian completion is trusted. The A map must not be reconstructed from \(A^\dagger A\); the cross Gram \(A^\dagger B\) requires the raw complex maps.

### Exact algebra already fixed, physical map still open

With the current convention
\[
K_{\rm comm}=[V,H_E^{\rm sine}],\qquad K_{\rm comm}^\dagger=-K_{\rm comm},
\]
the two-leg conversion is
\[
B_H=-\hbar^{-2}B_{\rm comm}.
\]
For the symmetric completion,
\[
B_{H,\rm sym}=-\frac{1}{2\hbar^2}
(B_{\rm comm}+B_{\rm comm}^\dagger).
\]
This is exact algebra conditional on the canonical commutator convention. It does not identify all factors in the physical extrinsic-curvature operator.

The existing classical regression fixes only the derivative-free relative kinetic-sector identity
\[
H_E^{\rm kin}=-\beta^2Q_{DW},\quad
H_L^{\rm corr}=(1+\beta^2)Q_{DW},\quad
H_E^{\rm kin}+H_L^{\rm corr}=Q_{DW}.
\]
It does not fix the absolute W5 normalization. The latter must be derived from a frozen canonical and regulator convention, including \(\kappa,\beta,\hbar,\ell_P\), regulator length, representation trace, orientation, lapse, and inverse-volume factors, with no double counting of factors already inside the computed map.

## Task C — joint kernel and physical causal cone

The joint-Gram gate is a finite W5 linear-algebra tool. Its output can establish only the numerical kernel dimension of the declared finite maps under the declared tolerance. It cannot turn a raw KKV candidate into a physical Hamiltonian.

The gate was hardened in this pass:
- all three output-basis identifiers must be explicit and equal;
- \(c_E,c_L\) and a normalization provenance string are mandatory;
- non-finite coefficients and empty output carriers are rejected;
- the five-dimensional right-kernel count includes structural zeros for rectangular output maps;
- regression tests cover rectangular nullity, exact test-map cancellation, implicit-normalization rejection, output-basis mismatch, empty carriers, and non-finite coefficients.

A physical joint kernel remains blocked until the common carrier and absolute coefficient convention are frozen. The physical causal cone is a still stronger result. It additionally requires a theory-derived clock or boundary-history ordering, a physical rigging map/inner product, connected physical response, a retarded Green function, and a universal IR comparison of scalar, vector, tensor, and matter sectors against an independently reconstructed spectral metric. The finite C8 relational positive control and the reduced free-fermion cone do not satisfy these requirements.

## Prior evidence that must not be overwritten

The saved V04 cross-branch W5 artifact contains the five-node Euclidean Gram and explicitly says: Euclidean node \(H_{E,v}\) only; no full Lorentzian H, D, or clock. It is a useful Euclidean anchor, not a Lorentzian amplitude source.

The historical \`docs/numerical_continuation/reverse_adjoint_middle_inventory.json\` reports 96 prepared middle inputs, 856 distinct source keys, zero evaluated middle actions, incomplete full reverse \(H_L\) column, and uncertified full error. It records real partial preparation, not completed Lorentzian amplitudes.

## Fail-closed status

| Deliverable | Status |
|---|---|
| Existing V09/V10 computational engines identified | audited |
| V10 10/12/14 full-node runs | queued; no result yet |
| Raw all-node W5 \(A\) and \(B\) common carrier | not yet assembled |
| Direct-adjoint closure and cutoff convergence | not yet numerically established |
| Minimal symmetric completion | algebraically defined; production matrix pending |
| Unique physical Hermitian ordering | open |
| Absolute physical normalization | open |
| Physical \(H_E+H_L\) joint kernel | open |
| Physical retarded response and causal cone | open |

No result is promoted from a finite reference, partial artifact, or queued workflow to a physical BQG claim.
