#!/usr/bin/env python3
"""Exact middle-prefix discriminator for the preregistered Lorentzian epsilon operator.

For one ordered pair (b,c) of distinct source-node incident edges this gate
computes, for every auxiliary path (i,j,k) in {0,1}^3,

    C_b(K) C_c(V) |psi>

using exactly the same sine-Hermitian ordering, zero-aware volume patch,
Peter-Weyl cutoff and covariant C(K)/C(V) implementations as the frozen
24-term logical-return calculation.

The result is independent of the outer edge a.  Therefore if every one of the
eight middle states is exactly empty within the frozen pruning convention, both
ordered triples with this same (b,c) prefix vanish before the outer C_a(K) is
applied.  A nonzero prefix is serialized path-by-path so a later outer-action
job can continue from it without recomputing C(V) and the middle C(K).

Zero/nonzero is a scientific result, not a PASS target. PASS checks only the
same leakage/cutoff integrity conditions as the preregistered full calculation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def encode_state(state):
    rows=[]
    for key,amp in sorted(state.items(),key=lambda kv:repr(kv[0])):
        spins,Kother,J2,M2,K12,K34=key
        rows.append({
            "spins":list(spins),"Kother":list(Kother),"J2":int(J2),"M2":int(M2),
            "K12":int(K12),"K34":int(K34),"amp":[float(amp.real),float(amp.imag)]
        })
    return rows


def ordered_pairs(source_v):
    n=tuple(FULL.RAW.PW.NEIG[source_v])
    return tuple((b,c) for b in n for c in n if b!=c)


def run(source_v=0,input_index=0,b=None,c=None,pair_index=None):
    ZVM.patch_and_clear()
    basis=FULL.RAW.PW.basis_full_jhalf()
    if len(basis)!=32:
        raise RuntimeError(f"expected 32 logical inputs, found {len(basis)}")
    if not (0<=input_index<len(basis)):
        raise ValueError("input_index outside logical basis")

    pairs=ordered_pairs(source_v)
    if pair_index is not None:
        if not (0<=pair_index<len(pairs)):
            raise ValueError("pair_index outside 0..11")
        pb,pc=pairs[pair_index]
        if b is not None and int(b)!=pb: raise ValueError("b conflicts with pair_index")
        if c is not None and int(c)!=pc: raise ValueError("c conflicts with pair_index")
        b,c=pb,pc
    if b is None or c is None:
        raise ValueError("provide --pair-index or both --b and --c")
    b=int(b);c=int(c)
    if (b,c) not in pairs:
        raise ValueError(f"({b},{c}) is not an ordered pair of distinct source-node neighbors")

    initial=basis[input_index]
    psi=FULL.RAW.CV.gauss_to_covariant({initial:1+0j},source_v)
    global_diag={
        "CV_complete_basis_leakage":0.0,
        "CK_outer_complete_basis_leakage":0.0,
        "CK_internal_volume_sector_leakage":0.0,
        "CK_complete_charge_basis_leakage":0.0,
    }
    paths=[]; max_spin=0.0

    old,caches=FULL.install_sine_ordering()
    try:
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    s1,leakV=FULL.RAW.COMP.C_volume_component(psi,source_v,c,k,i,FULL.JMAX2)
                    global_diag["CV_complete_basis_leakage"]=max(global_diag["CV_complete_basis_leakage"],float(leakV))
                    if s1:
                        s2,d2=FULL.RAW.KCOMP.C_K_component(s1,source_v,b,j,k,FULL.JMAX2)
                        FULL.update_diag(global_diag,d2)
                    else:
                        s2={}; d2={}
                    max_spin=max(max_spin,FULL.max_spin(s1),FULL.max_spin(s2))
                    paths.append({
                        "indices":[i,j,k],
                        "after_CV_support":len(s1),
                        "after_CV_norm":FULL.norm(s1),
                        "after_middle_CK_support":len(s2),
                        "after_middle_CK_norm":FULL.norm(s2),
                        "middle_state":encode_state(s2),
                    })
        cache_info={name:{"hits":fun.cache_info().hits,"misses":fun.cache_info().misses,"currsize":fun.cache_info().currsize} for name,fun in caches.items()}
    finally:
        FULL.restore_ordering(old)

    zero=all(row["after_middle_CK_support"]==0 for row in paths)
    hard={
        "CV_complete_basis_leakage_below_1e-9":global_diag["CV_complete_basis_leakage"]<1e-9,
        "CK_outer_complete_basis_leakage_below_1e-9":global_diag["CK_outer_complete_basis_leakage"]<1e-9,
        "CK_internal_volume_sector_leakage_below_1e-9":global_diag["CK_internal_volume_sector_leakage"]<1e-9,
        "spin_cutoff_respected":max_spin<=FULL.JMAX2/2+1e-12,
        "exactly_eight_aux_paths":len(paths)==8,
    }
    return {
        "status":"exact Lorentzian middle-prefix discriminator",
        "passed":bool(all(hard.values())),
        "science_status":"MIDDLE_PREFIX_ZERO_PATHWISE" if zero else "MIDDLE_PREFIX_NONZERO",
        "source_node":source_v,"input_logical_basis_index":input_index,"input_K_labels":list(initial[1]),
        "Jmax":FULL.JMAX2/2,"ordered_pair":{"b":b,"c":c,"pair_index":pairs.index((b,c))},
        "definition":"all eight auxiliary states C_b(K) C_c(V)|psi> before the outer C_a(K)",
        "prefix_zero_pathwise":zero,
        "max_spin_reached":max_spin,
        "max_diagnostics":global_diag,
        "hard_integrity_checks":hard,
        "paths":paths,
        "runtime_exact_cache":cache_info,
        "implication":(
            "Both full ordered triples sharing this (b,c) prefix vanish before the outer C_a(K)."
            if zero else
            "At least one middle path is nonzero; only the corresponding outer C_a(K) continuations can decide the full ordered triples."
        ),
        "claim_boundary":"Execution/selection-rule diagnostic for one preregistered logical input only; no physical projector, pole, dark matter or dark energy claim follows."
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-node',type=int,default=0)
    ap.add_argument('--input-index',type=int,default=0)
    ap.add_argument('--pair-index',type=int)
    ap.add_argument('--b',type=int);ap.add_argument('--c',type=int)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();out=run(a.source_node,a.input_index,a.b,a.c,a.pair_index)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='paths'},indent=2))
    return 0 if out['passed'] else 1

if __name__=='__main__': raise SystemExit(main())
