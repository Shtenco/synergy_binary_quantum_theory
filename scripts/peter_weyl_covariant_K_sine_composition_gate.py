#!/usr/bin/env python3
"""Physical sine-ordered adapter for the generalized C_e(K) composition gate.

Reuses the already-enforced exact cached composition engine, but replaces H_E
consistently in BOTH independent representations by

    H_E^sine = (T - T^dagger)/(2 i).

Only operator ordering changes. Holonomy paths, zero-aware volume, projections,
cutoffs and all acceptance thresholds are inherited unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import peter_weyl_covariant_K_composition_cached_gate as BASE
import peter_weyl_covariant_K_composition_gate as KC
import peter_weyl_euclidean_sine_ordering_gate as SINE


def projection_charge_nodes(charged_nodes, spec):
    """Include all vertices carrying open charge in this graph-changing term.

    A T-sequence touches the four vertices in its oriented spec. Its open
    fundamental indices are not confined to the source vertex: projecting
    the other touched vertices to Gauss intertwiners discards genuine branch
    support. Preserve caller-declared charges and add the term-local support.
    """
    v, a, b, c = spec
    return tuple(dict.fromkeys(tuple(charged_nodes) + (v, a, b, c)))


def complete_HE_sine(key, source_v, Jmax2, charged_nodes=(0,1), projection_diagnostics=None):
    base=KC.CK.branch_from_key(key)
    out={}; max_v=0.0; max_b=0.0
    for sign,spec in KC.PW.oriented_specs(source_v):
        v,a,b,c=spec
        term_charged_nodes = projection_charge_nodes(charged_nodes, spec)
        for adj in (False,True):
            pref=(+0.5j if adj else -0.5j)*sign
            for coef,seq0 in KC.PW.T_sequences(v,a,b,c):
                seq=KC.PW.adjoint_sequence(seq0) if adj else seq0
                branches,vleak=KC.CK.apply_sequence_to_branch(base,seq,source_v,Jmax2)
                max_v=max(max_v,float(vleak))
                for br in branches:
                    projected,bleak=KC.CK.project_branch_complete_charges(
                        br, term_charged_nodes, diagnostics=projection_diagnostics
                    )
                    max_b=max(max_b,float(bleak))
                    KC.CK.add_state(out,projected,pref*coef)
    return out,max_v,max_b


def complete_HE_sine_state(state, source_v, Jmax2, charged_nodes=(0,1), projection_diagnostics=None):
    """Apply the actual branch/recoupling path to a complex sparse superposition.

    Unlike a post-hoc matrix multiplication, each input coefficient is placed
    on the initial branch before the graph-changing sequences and projection
    run. Coherent contributions from all inputs and oriented terms are then
    accumulated in the output sparse map.
    """
    if not state:
        raise ValueError("input superposition must be non-empty")
    out = {}
    max_v = 0.0
    max_b = 0.0
    for key, amplitude in state.items():
        amp = complex(amplitude)
        if not (math.isfinite(amp.real) and math.isfinite(amp.imag)):
            raise ValueError("input superposition contains a non-finite amplitude")
        if amp == 0j:
            continue
        base = KC.CK.branch_from_key(key, amp=amp)
        for sign, spec in KC.PW.oriented_specs(source_v):
            v, a, b, c = spec
            term_charged_nodes = projection_charge_nodes(charged_nodes, spec)
            for adj in (False, True):
                pref = (+0.5j if adj else -0.5j) * sign
                for coef, seq0 in KC.PW.T_sequences(v, a, b, c):
                    seq = KC.PW.adjoint_sequence(seq0) if adj else seq0
                    branches, vleak = KC.CK.apply_sequence_to_branch(
                        base, seq, source_v, Jmax2
                    )
                    max_v = max(max_v, float(vleak))
                    for branch in branches:
                        projected, bleak = KC.CK.project_branch_complete_charges(
                            branch, term_charged_nodes, diagnostics=projection_diagnostics
                        )
                        max_b = max(max_b, float(bleak))
                        KC.CK.add_state(out, projected, pref * coef)
    return out, max_v, max_b


def gauss_HE_sine_with_historical_K_cutoff(state,node,Jmax2):
    """Match peter_weyl_lorentzian_K_block_gate.apply_HE_local exactly.

    The historical Gauss K reference prunes H_E at 1e-9.  The standalone sine
    ordering equivalence audit keeps 1e-10 to expose representation tails, but
    C(K_sine) must inherit the same 1e-9 Gauss K cutoff as C(K_plus).
    """
    return KC.PW.prune_state(SINE.safe_H_sine(state,node,Jmax2),1e-9)


def run(v=0,w=1):
    old_complete=KC.CK.apply_HE_complete_key
    old_gauss=KC.KG.apply_HE_local
    KC.CK.apply_HE_complete_key=complete_HE_sine
    KC.KG.apply_HE_local=gauss_HE_sine_with_historical_K_cutoff
    if hasattr(KC.CK.HE_complete_cached,'cache_clear'):
        KC.CK.HE_complete_cached.cache_clear()
    try:
        out=BASE.run(v,w)
    finally:
        KC.CK.apply_HE_complete_key=old_complete
        KC.KG.apply_HE_local=old_gauss
        if hasattr(KC.CK.HE_complete_cached,'cache_clear'):
            KC.CK.HE_complete_cached.cache_clear()
    out['euclidean_ordering']='sine-Hermitian (T-T^dagger)/(2i)'
    out['K_definition']='K_sine=[V,H_E^sine]'
    out['gauss_HE_cutoff']=1e-9
    out['historical_Gauss_K_cutoff_preserved']=True
    out['old_plus_K_reused']=False
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--v',type=int,default=0)
    ap.add_argument('--w',type=int,default=1)
    ap.add_argument('--output',type=Path)
    a=ap.parse_args(); out=run(a.v,a.w); text=json.dumps(out,indent=2); print(text)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text+'\n',encoding='utf-8')
    return 0 if out.get('passed',False) else 1

if __name__=='__main__':
    raise SystemExit(main())
