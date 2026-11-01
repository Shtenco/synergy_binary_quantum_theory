#!/usr/bin/env python3
"""Noninvasive REPLAY of actual node1 C_b(K) saved source + leakage diagnostics.
Monkeypatch OBSERVATION ONLY on SU2 primitive-basis projection; compare entire
replayed 120-key output with the archived KB packet before asserting any fact.
Does not test physical full Lorentzian HDA or physical TT.
"""
import argparse,json,math,collections,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import bqg_lorentzian_checkpointed_triple as Base
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_covariant_K_composition_gate as KC
import peter_weyl_lorentzian_ordered_triple_gate as RAW

def run(root):
    channel=(0,0,0)
    # All source tags and baseline hashes are checked by Base.load_result.
    kb_meta=json.loads((root/"kb_000.json").read_text(encoding="utf-8"))
    fr=kb_meta["freeze"]
    for name,val in fr["modules"].items():
        if Base.digest(Path(__file__).resolve().parent/name)!=val:
            raise RuntimeError("UPSTREAM_CODE_CHANGED: "+name)
    cv,cvmeta=Base.load_result(root,"cv",channel,fr)
    saved,savedmeta=Base.load_result(root,"kb",channel,fr)
    counts=collections.Counter();mass=collections.Counter();noteworthy=[]
    LP.JMAX2=Base.JMAX2
    restore,cache=LP.install_sine_cached_stack()
    original=KC.CK.project_branch_complete_charges
    def observe(branch,charged_nodes=(0,1),tol=1e-11):
        out,leak=original(branch,charged_nodes,tol)
        amp=branch[2]
        weight=abs(amp)**2
        counts["total"]+=1;mass["raw_weight_sum"]+=weight
        if leak>1e-8:
            counts["nontrivial_leak"]+=1;mass["branch_weight_at_nonzero_leak"]+=weight
        if leak>.9:
            counts["near_total_leak"]+=1;mass["branch_weight_at_near_total_leak"]+=weight
            if len(noteworthy)<8:
                noteworthy.append(dict(leak=float(leak),amplitude_abs=float(abs(amp)),
                    source_spin2=list(branch[0]),charged_nodes=list(charged_nodes),
                    projected_basis_support=len(out)))
        return out,leak
    KC.CK.project_branch_complete_charges=observe
    try:
        computed,diag=RAW.KCOMP.C_K_component(cv,fr["source"],fr["triple"][1],0,0,fr["JMAX2"])
    finally:
        KC.CK.project_branch_complete_charges=original
        restore()
    discarded=len(computed)-sum(k[2] in (0,2) for k in computed)
    pruned={k:v for k,v in computed.items() if k[2] in (0,2)}
    d2=Base.norm({k:pruned.get(k,0)-saved.get(k,0) for k in set(pruned)|set(saved)})
    tol=1e-9*max(1,Base.norm(saved))
    result=dict(replay_status="PASS" if d2<tol else "FAIL",
        exact_input_source="archived CV stage",archived_KB_sha256=savedmeta["sha256"],
        archived_norm=Base.norm(saved),replayed_norm=Base.norm(pruned),archived_support=len(saved),
        replay_support=len(pruned),discarded_J2_count=discarded,
        absolute_replay_error=d2,relative_replay_error=d2/max(Base.norm(saved),1e-30),
        unprojected_primitive_leak_count=dict(counts),
        unprojected_primitive_weight=dict(mass),
        suspicious_branches=noteworthy,
        upstream_reported_CK_diagnostics=diag,
        statement="Primitive-branch leaks do NOT alone measure post-H_E-summed gauge/charge leakage; no physical gate promoted.")
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--checkpoint-dir",type=Path,required=True)
    p.add_argument("--output",type=Path); a=p.parse_args()
    result=run(a.checkpoint_dir)
    s=json.dumps(result,ensure_ascii=False,indent=2)
    print(s,flush=True)
    if a.output:a.output.write_text(s+"\n",encoding="utf-8")
    if result["replay_status"]!="PASS":raise SystemExit(2)
