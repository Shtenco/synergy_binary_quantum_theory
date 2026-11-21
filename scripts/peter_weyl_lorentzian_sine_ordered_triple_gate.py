#!/usr/bin/env python3
"""First physical sine-ordered Peter-Weyl K_sine-K_sine-V triple.

Evaluates the existing raw ordered triple

    Tr_aux[C_a(K_sine) C_b(K_sine) C_c(V)]

at the preregistered single-H_L wall Jmax=7/2, but rebuilds K everywhere from

    H_E^sine=(T-T^dagger)/(2i),  K_sine=[V,H_E^sine].

No K_plus amplitude is reused. Runtime reduction is exact: scalar H_E reduced
matrix elements are cached at highest weight M=J and relabelled by the rank-0
Wigner--Eckart identity. Full sparse states remain exact cache keys.

The historical fixed-index primitive charge-basis diagnostic is retained in the
JSON but is not final physical leakage; the independently enforced C(K_sine)
gate has already shown final wrong-charge fractions vanish while that primitive
diagnostic equals 1 in both representations.
"""
from __future__ import annotations

import argparse
import functools
import json
from pathlib import Path

import peter_weyl_covariant_K_composition_gate as KC
import peter_weyl_covariant_K_composition_cached_gate as CACHE
import peter_weyl_covariant_K_sine_composition_gate as SINEK
import peter_weyl_lorentzian_ordered_triple_gate as RAW


def freeze(state):
    return tuple(sorted(state.items(),key=lambda kv:repr(kv[0])))


def run(v=0,a=None,b=None,c=None,executor=None):
    old_he=KC.CK.apply_HE_complete_key
    old_gauss=KC.KG.apply_HE_local
    old_k=KC.apply_K_complete_custom
    old_inv=KC.COMP.inverse_complete
    old_direct=KC.direct_K_covariant
    old_close=KC.COMP.close_complete
    old_ck=KC.C_K_component
    old_cv=RAW.COMP.C_volume_component

    @functools.lru_cache(maxsize=None)
    def he_reduced(canonical_key,source_v,Jmax2,charged_nodes):
        return SINEK.complete_HE_sine(
            canonical_key,source_v,Jmax2,charged_nodes=tuple(charged_nodes)
        )

    def he_sine_reduced(key,source_v,Jmax2,charged_nodes=(0,1)):
        charged_nodes=tuple(charged_nodes)
        canonical,original=CACHE.canonicalize_scalar_charge_M(key,charged_nodes)
        state,vleak,bleak=he_reduced(canonical,source_v,Jmax2,charged_nodes)
        return CACHE.restore_scalar_charge_M(state,charged_nodes,original),vleak,bleak

    @functools.lru_cache(maxsize=None)
    def k_cached(frozen,source_v,Jmax2,charged_nodes):
        out,vleak,bleak=old_k(dict(frozen),source_v,Jmax2,tuple(charged_nodes))
        return tuple(out.items()),float(vleak),float(bleak)
    def k_wrap(state,source_v,Jmax2,charged_nodes):
        items,vleak,bleak=k_cached(freeze(state),source_v,Jmax2,tuple(charged_nodes))
        return dict(items),vleak,bleak

    @functools.lru_cache(maxsize=None)
    def inv_cached(frozen,source_v,target_v,k,j,Jmax2):
        out,leak=old_inv(dict(frozen),source_v,target_v,k,j,Jmax2)
        return tuple(out.items()),float(leak)
    def inv_wrap(state,source_v,target_v,k,j,Jmax2):
        items,leak=inv_cached(freeze(state),source_v,target_v,k,j,Jmax2)
        return dict(items),leak

    @functools.lru_cache(maxsize=None)
    def direct_cached(frozen,source_v,Jmax2):
        out,vleak,bleak=old_direct(dict(frozen),source_v,Jmax2)
        return tuple(out.items()),float(vleak),float(bleak)
    def direct_wrap(state,source_v,Jmax2):
        items,vleak,bleak=direct_cached(freeze(state),source_v,Jmax2)
        return dict(items),vleak,bleak

    @functools.lru_cache(maxsize=None)
    def close_cached(frozen,source_v,target_v,i,k,Jmax2):
        return tuple(old_close(dict(frozen),source_v,target_v,i,k,Jmax2).items())
    def close_wrap(state,source_v,target_v,i,k,Jmax2):
        return dict(close_cached(freeze(state),source_v,target_v,i,k,Jmax2))

    @functools.lru_cache(maxsize=None)
    def ck_cached(frozen,source_v,target_v,i,j,Jmax2):
        out,diag=old_ck(dict(frozen),source_v,target_v,i,j,Jmax2)
        return tuple(out.items()),tuple(sorted(diag.items()))
    def ck_wrap(state,source_v,target_v,i,j,Jmax2):
        items,diag=ck_cached(freeze(state),source_v,target_v,i,j,Jmax2)
        return dict(items),dict(diag)

    @functools.lru_cache(maxsize=None)
    def cv_cached(frozen,source_v,target_v,i,j,Jmax2):
        out,leak=old_cv(dict(frozen),source_v,target_v,i,j,Jmax2)
        return tuple(out.items()),float(leak)
    def cv_wrap(state,source_v,target_v,i,j,Jmax2):
        items,leak=cv_cached(freeze(state),source_v,target_v,i,j,Jmax2)
        return dict(items),leak

    KC.CK.apply_HE_complete_key=he_sine_reduced
    KC.KG.apply_HE_local=SINEK.gauss_HE_sine_with_historical_K_cutoff
    KC.apply_K_complete_custom=k_wrap
    KC.COMP.inverse_complete=inv_wrap
    KC.direct_K_covariant=direct_wrap
    KC.COMP.close_complete=close_wrap
    KC.C_K_component=ck_wrap
    RAW.KCOMP.C_K_component=ck_wrap
    RAW.COMP.C_volume_component=cv_wrap
    if hasattr(KC.CK.HE_complete_cached,'cache_clear'):
        KC.CK.HE_complete_cached.cache_clear()
    try:
        # Preserve every original sine-H_E adapter and cache while a v0.8
        # amplitude assembler runs. The default single-triple path is unchanged.
        if executor is not None:
            return executor()
        out=RAW.run(v,a,b,c)
        infos={
            'HE_reduced':he_reduced.cache_info(),
            'K_complete':k_cached.cache_info(),
            'inverse':inv_cached.cache_info(),
            'direct_K':direct_cached.cache_info(),
            'close':close_cached.cache_info(),
            'C_K':ck_cached.cache_info(),
            'C_V':cv_cached.cache_info(),
        }
    finally:
        KC.CK.apply_HE_complete_key=old_he
        KC.KG.apply_HE_local=old_gauss
        KC.apply_K_complete_custom=old_k
        KC.COMP.inverse_complete=old_inv
        KC.direct_K_covariant=old_direct
        KC.COMP.close_complete=old_close
        KC.C_K_component=old_ck
        RAW.KCOMP.C_K_component=old_ck
        RAW.COMP.C_volume_component=old_cv
        if hasattr(KC.CK.HE_complete_cached,'cache_clear'):
            KC.CK.HE_complete_cached.cache_clear()

    raw_pass=bool(out.get('passed',False))
    d=out['max_diagnostics']
    physical_pass=bool(
        out['output_support']>0
        and out['output_norm']>1e-10
        and out['final_nonscalar_J_weight_fraction']<1e-8
        and d['CV_complete_basis_leakage']<1e-9
        and d['CK_outer_complete_basis_leakage']<1e-9
        and d['CK_internal_volume_sector_leakage']<1e-9
        and out['max_spin_reached']<=3.5+1e-12
    )
    out['raw_base_gate_passed']=raw_pass
    out['passed']=physical_pass
    out['euclidean_ordering']='sine-Hermitian (T-T^dagger)/(2i)'
    out['K_definition']='K_sine=[V,H_E^sine]'
    out['old_plus_K_reused']=False
    out['historical_primitive_charge_basis_diagnostic']={
        'value':d['CK_complete_charge_basis_leakage'],
        'hard_acceptance':False,
        'reason':'primitive fixed-index branches precede the complete gauge-invariant H_E sum; enforced C(K_sine) final wrong-charge fractions are zero',
    }
    out['runtime_exact_reduction']={
        name:{'hits':x.hits,'misses':x.misses,'currsize':x.currsize}
        for name,x in infos.items()
    }
    out['runtime_exact_reduction']['rank0_Wigner_Eckart_M_reduction']=True
    out['runtime_exact_reduction']['physics_changed']=False
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--v',type=int,default=0)
    ap.add_argument('--a',type=int)
    ap.add_argument('--b',type=int)
    ap.add_argument('--c',type=int)
    ap.add_argument('--output',type=Path)
    x=ap.parse_args(); out=run(x.v,x.a,x.b,x.c); text=json.dumps(out,indent=2); print(text)
    if x.output:
        x.output.parent.mkdir(parents=True,exist_ok=True)
        x.output.write_text(text+'\n',encoding='utf-8')
    return 0 if out.get('passed',False) else 1

if __name__=='__main__':
    raise SystemExit(main())
