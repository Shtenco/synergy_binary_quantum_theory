# Common Gauss-key habitat audit: Euclidean and Lorentzian BQG maps

**Status:** bounded finite operator-domain compatibility experiment. This is a prerequisite check for BH-0, not a physical projector result.

Reproduction:

```bash
python scripts/bqg_common_gauss_habitat_audit.py
```

The CI default uses the five frozen W5 seed labels (K\in\{0,2,4,6,8\}), source node (v=0), and (J_{max2}=10).

## 1. Actual operators applied to the same seed

For each W5 Gauss-key seed

[
|s_K\rangle=|((4)^{10},(K,0,0,0,0))\rangle,
]

the script applies the existing Euclidean graph-changing map

[
E_v|s_K\rangle=H_v^E|s_K\rangle
]

from `k5_peter_weyl_safe_hda_column.py), and the existing raw Lorentzian KKV map plus its independently ordered direct adjoint,

[
L_v|s_K\rangle,\qquad L_v^{\dagger,\mathrm{direct}}|s_K\rangle,
]

using the exact finite covariant-volume and sine-(K) components already used by the W5 closure runner. It uses the same zero-aware volume-migration convention as that runner.

The Lorentzian outputs are converted back to Gauss scalar keys with the existing exact scalar map. No arbitrary dense matrix or hand-selected common kernel is introduced.

## 2. Common kinematical key space

All three outputs are represented in the same key schema

[
(\mathbf j,\mathbf K),\qquad
\mathbf j\in\mathbb Z_{\ge0}^{10},\quad
\mathbf K\in\mathbb Z_{\ge0}^{5}.
]

For each seed the audit records:
- output support sizes and norms;
- intersection and union of Euclidean/Lorentzian supports;
- normalized Euclidean–Lorentzian output overlap;
- (j=0) cylindrical graph-change norm fractions;
- direct-adjoint output support;
- non-scalar covariant output rejected by the Gauss-scalar map;
- component leakage diagnostics supplied by the underlying operators.

The normalized overlap is a diagnostic of shared output support, not an assertion that the two operators are equal or that their commutator closes.

## 3. Exact pass boundary

A PASS means only:
1. the frozen W5 seeds and all retained outputs have the same Gauss-key data schema;
2. the three actual finite maps return nonzero columns on the tested seeds;
3. the declared Lorentzian scalar reduction rejects no residual output above tolerance.

A PASS does **not** prove:
- that (H_E) and raw (L) are members of one anomaly-free constraint family;
- that the direct Lorentzian adjoint is the chosen physical Hermitian ordering;
- that the combined operator domain is closed under every constraint and adjoint;
- that a common kernel exists or is nontrivial;
- that a finite master constraint made from these maps is the physical BQG projector;
- positivity of the physical inner product or existence of a refinement/rigging-map limit.

In particular, the finite Euclidean three-node HDA control, the W5 Lorentzian adjoint-closure control and the C8 relational positive control remain distinct constructions. This audit checks a shared *kinematical seed/key interface*; it does not identify their dynamical labels or clocks.

## 4. Next hard gate

If this audit passes, the next step is a bounded joint-domain closure experiment:

1. Freeze one common graph, node set, spin cutoff, domains, ordering and inner product.
2. Construct actual matrix columns for every declared constraint and adjoint on the same basis.
3. Expand the basis by every nonzero output until closure, or fail closed with an explicit omitted-support certificate.
4. Verify the declared constraint algebra on that enlarged domain.
5. Only then form (M_G=\sum_{A,B}C_A^\dagger G^{AB}C_B), compare its common kernel across positive (G), and study a second regulator.

The current script deliberately stops before step 2 becomes a joint matrix construction. No physical projector is claimed.
