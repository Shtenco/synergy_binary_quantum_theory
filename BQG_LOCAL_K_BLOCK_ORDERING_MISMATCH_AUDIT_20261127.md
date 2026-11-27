# Local K-block implementation audit — 2026-11-27

## Finding: the existing local K script is not the sine-Hermitian K block

The normalization note `BQG_W5_K_COMMUTATOR_NORMALIZATION.md` defines
[
H_E^{sine}=(T-T^\dagger)/(2i),\qquad K_{comm}=[V,H_E^{sine}].
]
However, `scripts/peter_weyl_lorentzian_K_block_gate.py` implements its Euclidean leg through:

```python
def apply_HE_local(state, v, Jmax2):
    return PW.prune_state(PW.apply_H_cached_state(state, v, Jmax2), 1e-9)
```

The imported `PW.apply_H_cached_state` in `scripts/k5_peter_weyl_safe_hda_column.py` constructs
[
H_E^{plus}=\tfrac12(T+T^\dagger)
]
by adding the raw T map and its adjoint with coefficient (+1/2) each (lines 284–290). It does **not** construct the sine combination ((T-T^\dagger)/(2i)).

Therefore the output of `peter_weyl_lorentzian_K_block_gate.py` is a commutator ([V,H_E^{plus}]), not the (K_{comm}) whose sign/(\hbar) normalization is derived in `BQG_W5_K_COMMUTATOR_NORMALIZATION.md`. The sine-normalization formula is algebraically correct, but it cannot be attached to the local K-block output until the implementation uses the same sine operator.

This is a code/documentation mismatch, not a claim that either Euclidean ordering is physically preferred. No coefficient or sign has been fitted.

## Additional numerical-scope limitations

The local K script also:
- prunes (V|psi\rangle) at (10^{-10});
- prunes (H_E|psi\rangle) at (10^{-9});
- prunes the commutator output at (10^{-9});
- uses lower-level cutoffs in the imported Peter–Weyl engine (including (10^{-12}), (10^{-11}), and (10^{-13}));
- checks reverse matrix elements only for the largest `reverse_samples` outputs (default 8), rather than assembling and testing the full finite matrix.

Thus the label “exact” in the script title/status should be read as “built from the exact finite recoupling/volume definitions within the declared implementation,” not as an untruncated exact operator identity or a full anti-Hermiticity certificate. The reported reverse-sample test is a finite spot check.

## Required repair sequence (audit branch only)

1. Add an explicit Euclidean-ordering selector to the local K gate, or create a separate sine-specific runner. The sine route must call the already existing sine-Hermitian complete Euclidean engine, not the (H_E^{plus}) helper.
2. Freeze one input state, cutoff, node, basis convention, and all internal cutoffs. Record the full cutoff ledger.
3. Build both ([V,H_E^{sine}]) and its independently ordered reverse matrix elements on the **same** finite reachable basis. Do not infer anti-Hermiticity from eight largest samples.
4. Separately build ([V,H_E^{plus}]) as a control. Do not mix its amplitudes with the sine result.
5. Check the exact algebraic relation (K_H=(i\hbar)^{-1}K_{comm}) only for the sine-defined anti-Hermitian commutator after its matrix adjoint has been verified.
6. For two Hermitian K legs, preserve the factor (-\hbar^{-2}). This remains only the conversion of the two chosen commutator legs; the full Thiemann/((1+\beta^2)), (\kappa), regulator and canonical normalization must be derived separately.

## Status

- Sine commutator normalization identity: **algebraically valid, conditional on the sine-defined anti-Hermitian input**.
- Existing `peter_weyl_lorentzian_K_block_gate.py` as a sine-K certificate: **NOT VALID; it currently calls the plus-Hermitian engine**.
- Existing local K amplitude block as a finite numerical prerequisite: **present but thresholded and only spot-checked for anti-Hermiticity**.
- Full (K_{ab}), physical Lorentzian normal generator, common-domain HDA and BQG null cone: **OPEN**.

This audit does not modify the production operator or claim a physical result.


## Precision correction to the normalization interpretation (same-day amendment)

The ordering mismatch above must be stated precisely:

- The sine-specific normalization document defines (K_{comm}^{sine}=[V,H_E^{sine}]).
- The existing local gate computes a different candidate (K_{comm}^{plus}=[V,H_E^{plus}]).
- **The algebraic conversion itself is not sine-specific.** For either choice, if (V) and the selected (H_E) are Hermitian on the same domain, then the exact commutator is anti-Hermitian and
  [
  K_H=(i\hbar)^{-1}K_{comm}
  ]
  is Hermitian. Consequently, the two-leg factor (-\hbar^{-2}) also applies to two such converted legs, whether their common Hermitian Euclidean input is (H_E^{sine}) or (H_E^{plus}).

The actual defect is therefore **operator-definition mismatch**, not invalidity of the universal anti-Hermitian-to-Hermitian conversion. The existing plus-ordered block may be normalized by the same algebraic identity after its full finite-matrix anti-Hermiticity is established. It must, however, be labelled (K^{plus}), and it cannot be reported as the sine-defined block in the sine normalization document.

Because the current gate prunes intermediate states and checks only eight reverse matrix elements, its numerical output has not yet established the exact finite-matrix premise. This amendment supersedes the overly strong sentence above that said the conversion “cannot be attached” to the local output. It can be attached conditionally as a plus-ordered conversion, with that explicit label and the same-domain adjoint check still required.
