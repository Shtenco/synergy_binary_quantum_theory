#!/usr/bin/env python3
"""Compute all 12 Lorentzian middle prefixes for one source/input in one cache session.

This is algebraically the same discriminator as
peter_weyl_lorentzian_middle_prefix_gate.py, but removes cross-prefix duplicate
work.  For fixed source/input,

    C_c(V)_{k i}|psi>

depends only on c,k,i, so there are 4*2*2 = 16 unique volume-leg states, not
one recomputation per (b,c,j).  The sine-ordering caches are also shared across
all 12 ordered prefixes.  Every emitted prefix JSON preserves the original
middle-prefix contract and can be consumed by the existing outer workers.

No zero/nonzero pattern is assumed in advance; all 12 prefixes are measured for
each new logical input.
"""
from __future__ import annotations

import argparse, json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_middle_prefix_gate as MID
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def run(source_v=0,input_index=0):
    ZVM.patch_and_clear()
    basis=FULL.RAW.PW.basis_full_jhalf()
    if len(basis)!=32: raise RuntimeError(f"expected 32 logical inputs, found {len(basis)}")
    if not (0<=input_index<len(basis)): raise ValueError("input_index outside logical basis")
    initial=basis[input_index]
    psi=FULL.RAW.CV.gauss_to_covariant({initial:1+0j},source_v)
    pairs=MID.ordered_pairs(source_v)
    neighbors=tuple(FULL.RAW.PW.NEIG[source_v])
    if len(neighbors)!=4 or len(pairs)!=12: raise RuntimeError("expected four-valent source / 12 ordered prefixes")

    # Unique C(V) states.  The original per-prefix loop repeats these across b
    # and across j, although C_c(V)_{ki} contains neither b nor j.
    cv_states={}; cv_leak={}; prefixes=[]
    old,caches=FULL.install_sine_ordering()
    try:
        for c in neighbors:
            for i in range(2):
                for k in range(2):
                    s1,leak=FULL.RAW.COMP.C_volume_component(psi,source_v,c,k,i,FULL.JMAX2)
                    cv_states[(c,i,k)]=s1
                    cv_leak[(c,i,k)]=float(leak)

        for pair_index,(b,c) in enumerate(pairs):
            diag={
                "CV_complete_basis_leakage":max(cv_leak[(c,i,k)] for i in range(2) for k in range(2)),
                "CK_outer_complete_basis_leakage":0.0,
                "CK_internal_volume_sector_leakage":0.0,
                "CK_complete_charge_basis_leakage":0.0,
            }
            paths=[]; max_spin=0.0
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        s1=cv_states[(c,i,k)]
                        if s1:
                            s2,d2=FULL.RAW.KCOMP.C_K_component(s1,source_v,b,j,k,FULL.JMAX2)
                            FULL.update_diag(diag,d2)
                        else:
                            s2={}
                        max_spin=max(max_spin,FULL.max_spin(s1),FULL.max_spin(s2))
                        paths.append({
                            "indices":[i,j,k],
                            "after_CV_support":len(s1),
                            "after_CV_norm":FULL.norm(s1),
                            "after_middle_CK_support":len(s2),
                            "after_middle_CK_norm":FULL.norm(s2),
                            "middle_state":MID.encode_state(s2),
                        })
            zero=all(r["after_middle_CK_support"]==0 for r in paths)
            hard={
                "CV_complete_basis_leakage_below_1e-9":diag["CV_complete_basis_leakage"]<1e-9,
                "CK_outer_complete_basis_leakage_below_1e-9":diag["CK_outer_complete_basis_leakage"]<1e-9,
                "CK_internal_volume_sector_leakage_below_1e-9":diag["CK_internal_volume_sector_leakage"]<1e-9,
                "spin_cutoff_respected":max_spin<=FULL.JMAX2/2+1e-12,
                "exactly_eight_aux_paths":len(paths)==8,
            }
            prefixes.append({
                "status":"exact Lorentzian middle-prefix discriminator",
                "execution_mode":"all_12_prefixes_shared_cache_v1",
                "passed":bool(all(hard.values())),
                "science_status":"MIDDLE_PREFIX_ZERO_PATHWISE" if zero else "MIDDLE_PREFIX_NONZERO",
                "source_node":int(source_v),
                "input_logical_basis_index":int(input_index),
                "input_K_labels":list(initial[1]),
                "Jmax":FULL.JMAX2/2,
                "ordered_pair":{"b":int(b),"c":int(c),"pair_index":int(pair_index)},
                "definition":"all eight auxiliary states C_b(K) C_c(V)|psi> before the outer C_a(K)",
                "prefix_zero_pathwise":bool(zero),
                "max_spin_reached":max_spin,
                "max_diagnostics":diag,
                "hard_integrity_checks":hard,
                "paths":paths,
                "implication":(
                    "Both full ordered triples sharing this (b,c) prefix vanish before the outer C_a(K)."
                    if zero else
                    "At least one middle path is nonzero; only the corresponding outer C_a(K) continuations can decide the full ordered triples."
                ),
                "claim_boundary":"Execution/selection-rule diagnostic for one preregistered logical input only; no physical projector or cosmological claim follows.",
            })
        cache_info={name:{"hits":fun.cache_info().hits,"misses":fun.cache_info().misses,"currsize":fun.cache_info().currsize} for name,fun in caches.items()}
    finally:
        FULL.restore_ordering(old)

    # Attach the common final cache audit after all prefixes have shared it.
    for p in prefixes: p["runtime_exact_cache"]=cache_info
    zero_indices=[p["ordered_pair"]["pair_index"] for p in prefixes if p["prefix_zero_pathwise"]]
    live_indices=[p["ordered_pair"]["pair_index"] for p in prefixes if not p["prefix_zero_pathwise"]]
    summary_checks={
        "twelve_prefixes_once":len(prefixes)==12 and {p["ordered_pair"]["pair_index"] for p in prefixes}==set(range(12)),
        "sixteen_unique_CV_states":len(cv_states)==16,
        "all_prefix_integrity_passed":all(p["passed"] for p in prefixes),
        "zero_live_partition_complete":len(zero_indices)+len(live_indices)==12 and not (set(zero_indices)&set(live_indices)),
    }
    return {
        "schema":"BQG_LORENTZIAN_ALL_MIDDLE_PREFIXES_V1",
        "passed":bool(all(summary_checks.values())),
        "source_node":int(source_v),"input_logical_basis_index":int(input_index),"input_K_labels":list(initial[1]),
        "Jmax":FULL.JMAX2/2,
        "unique_CV_state_count":len(cv_states),
        "naive_separate_prefix_CV_requests":96,
        "explicit_unique_CV_evaluations":16,
        "zero_prefix_indices":zero_indices,"nonzero_prefix_indices":live_indices,
        "zero_prefix_count":len(zero_indices),"nonzero_prefix_count":len(live_indices),
        "summary_checks":summary_checks,"runtime_exact_cache":cache_info,"prefixes":prefixes,
        "claim_boundary":"Batched exact middle-prefix packet for one source/input. It classifies which outer continuations are required but is not an H_L column, master, HDA certificate or physical projector."
    }


def write_outputs(out,output_dir:Path):
    output_dir.mkdir(parents=True,exist_ok=True)
    summary={k:v for k,v in out.items() if k!="prefixes"}
    (output_dir/'middle_prefix_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    for p in out['prefixes']:
        idx=p['ordered_pair']['pair_index']
        (output_dir/f'prefix_{idx}.json').write_text(json.dumps(p,indent=2)+'\n',encoding='utf-8')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-node',type=int,default=0);ap.add_argument('--input-index',type=int,default=0)
    ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    out=run(a.source_node,a.input_index);write_outputs(out,a.output_dir)
    print(json.dumps({k:v for k,v in out.items() if k!='prefixes'},indent=2))
    return 0 if out['passed'] else 1

if __name__=='__main__': raise SystemExit(main())
