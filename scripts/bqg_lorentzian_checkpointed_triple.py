#!/usr/bin/env python3
"""Resumable genuine full-state BQG Lorentzian triple, with 8 aux channels x 3 stages.

NO substitution of model amplitudes. Same LP/R AW/CV/KCOMP operator kernels,
Jmax=7/2, charged-node rules, scalar pruning and scalar closure as the original
peter_weyl_lorentzian_fullstate_triple_worker.py. Checkpoints include ALL
covariant keys (not just logical states) and max-leak diagnostics.
Stage options: cv, kb, ka, finish. Each stage is one durable GitHub Actions job.
This computes a *single raw ordered triple*, NOT the full Hermitian Lorentzian H.
"""
from __future__ import annotations
import argparse, ast, hashlib, itertools, json, math, os, sys, time
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_lorentzian_gauss_action_gate as LGA
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_lorentzian_ordered_triple_gate as RAW

JMAX2=7
TOL=1e-10
LEAK_FIELDS=("CV_complete_basis_leakage","CK_outer_complete_basis_leakage",
             "CK_internal_volume_sector_leakage","CK_complete_charge_basis_leakage")
STAGES=("cv","kb","ka","finish")
CHECKPOINT_SCHEMA = 2
ROOT_MODULES = (
    "bqg_lorentzian_checkpointed_triple",
    "bqg_lorentzian_ka_partition",
)

def runtime_fingerprint():
    return {"python": sys.version.split()[0], "numpy": np.__version__}

def _local_module_path(name):
    relative = Path(*name.split("."))
    module_file = HERE / relative.with_suffix(".py")
    package_file = HERE / relative / "__init__.py"
    if module_file.is_file():
        return module_file
    if package_file.is_file():
        return package_file
    return None

def source_hashes():
    """Hash the complete transitive closure of local imports used by both workers."""
    pending = list(ROOT_MODULES)
    visited = set()
    hashes = {}
    while pending:
        name = pending.pop()
        if name in visited:
            continue
        visited.add(name)
        path = _local_module_path(name)
        if path is None:
            continue
        relative = path.relative_to(HERE).as_posix()
        hashes[relative] = digest(path)
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=relative)
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                pending.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                if node.level == 0 and node.module:
                    pending.append(node.module)
                elif node.level > 0:
                    package = path.relative_to(HERE).with_suffix("").parts[:-1]
                    keep = max(0, len(package) - node.level + 1)
                    parts = package[:keep]
                    if node.module:
                        parts += tuple(node.module.split("."))
                    if parts:
                        pending.append(".".join(parts))
    return dict(sorted(hashes.items()))
def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        while True:
            block=f.read(1024*1024)
            if not block:break
            h.update(block)
    return h.hexdigest()
def write_atomic(path,write):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+".partial")
    with tmp.open("wb") as f:
        write(f);f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)
def persist_state(path,state):
    rows=sorted(state.items(),key=lambda t:repr(t[0]))
    n=len(rows)
    spins=np.array([key[0] for key,v in rows],dtype=np.int16).reshape(n,len(PW.EDGES))
    others=np.array([key[1] for key,v in rows],dtype=np.int16).reshape(n,len(PW.VERT))
    ints=np.array([key[2:] for key,v in rows],dtype=np.int16).reshape(n,4)
    amps=np.asarray([v for key,v in rows],dtype=np.complex128)
    write_atomic(path,lambda f:np.savez_compressed(f,spins=spins,Kother=others,
        charges=ints,amp=amps))
def restore_state(path):
    with np.load(path,allow_pickle=False) as z:
        spins=z["spins"];others=z["Kother"];ints=z["charges"];amplitudes=z["amp"]
    if not(len(spins)==len(others)==len(ints)==len(amplitudes)):raise ValueError("corrupt row alignment")
    if spins.shape[1:]!=(len(PW.EDGES),) or others.shape[1:]!=(len(PW.VERT),) or ints.shape[1:]!=(4,):
        raise ValueError("corrupt covariant key dimensions")
    out={}
    for s,k,j,v in zip(spins,others,ints,amplitudes):
        key=(tuple(map(int,s)),tuple(map(int,k)),*map(int,j))
        if key in out:raise ValueError("duplicate covariant key")
        if not (np.isfinite(v.real) and np.isfinite(v.imag)):raise ValueError("nonfinite amplitude")
        out[key]=complex(v)
    return out
def norm(state):
    return math.sqrt(sum(abs(z)**2 for z in state.values()))
def max_spin(state):
    return max((max(k[0]) for k in state),default=0)/2
def selected(args):
    all_=list(itertools.product(range(2),repeat=3))
    return all_ if args.channel=="all" else [tuple(map(int,args.channel))]
def prefix(d,stage,ijk):
    return d/(stage+"_"+''.join(map(str,ijk)))
def validate_leakage_map(diag, require_all=False, context="checkpoint"):
    if not isinstance(diag, dict):
        raise ValueError(f"{context}: leakage_max must be an object")
    unknown = set(diag) - set(LEAK_FIELDS)
    if unknown:
        raise ValueError(f"{context}: unknown leakage fields {sorted(unknown)}")
    missing = set(LEAK_FIELDS) - set(diag)
    if require_all and missing:
        raise ValueError(f"{context}: missing required leakage diagnostics {sorted(missing)}")
    clean = {}
    for key, value in diag.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{context}: {key} must be numeric")
        value = float(value)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{context}: {key} must be finite and nonnegative")
        clean[key] = value
    return clean

def validate_operator_provenance(meta, context="checkpoint"):
    if not isinstance(meta, dict):
        raise ValueError(f"{context}: operator provenance must be an object")
    input_support = meta.get("operator_input_support")
    if isinstance(input_support, bool) or not isinstance(input_support, int) or input_support < 0:
        raise ValueError(f"{context}: invalid operator input support")
    applied = meta.get("operator_stage_applied")
    if not isinstance(applied, bool) or applied is not (input_support > 0):
        raise ValueError(f"{context}: operator-stage execution flag mismatch")
    expected_fields = {
        "CK_outer_complete_basis_leakage",
        "CK_internal_volume_sector_leakage",
        "CK_complete_charge_basis_leakage",
    }
    expected = expected_fields if input_support == 0 else set()
    actual = meta.get("exact_empty_input_zero_fields")
    if not isinstance(actual, list) or set(actual) != expected or len(actual) != len(expected):
        raise ValueError(f"{context}: empty-input diagnostic provenance mismatch")
    if meta.get("diagnostic_provenance_scope") != "current_operator_stage":
        raise ValueError(f"{context}: diagnostic provenance scope is missing or invalid")
    stage_diag = validate_leakage_map(
        meta.get("stage_leakage_diagnostics"), require_all=True,
        context=f"{context} current operator stage",
    )
    if set(stage_diag) != expected_fields:
        raise ValueError(f"{context}: current stage must report all three CK leakage channels")
    if input_support == 0 and any(stage_diag[field] != 0.0 for field in expected_fields):
        raise ValueError(f"{context}: exact-empty-input stage diagnostics must be exactly zero")
    return {"operator_input_support": input_support,
            "operator_stage_applied": applied,
            "diagnostic_provenance_scope": "current_operator_stage",
            "exact_empty_input_zero_fields": sorted(expected),
            "stage_leakage_diagnostics": stage_diag}


def validate_checkpoint_chain(cv_meta, kb_meta, ka_meta, cv_support, kb_support):
    """Validate provenance and max aggregation across cv -> kb -> ka."""
    cv_diag = validate_leakage_map(
        cv_meta.get("leakage_max"), require_all=False, context="cv checkpoint"
    )
    if set(cv_diag) != {"CV_complete_basis_leakage"}:
        raise ValueError("cv checkpoint must contain exactly the CV leakage diagnostic")

    kb_provenance = validate_operator_provenance(kb_meta, context="kb checkpoint")
    ka_provenance = validate_operator_provenance(ka_meta, context="ka checkpoint")
    if kb_provenance["operator_input_support"] != cv_support:
        raise ValueError("kb operator input support does not match cv output support")
    if ka_provenance["operator_input_support"] != kb_support:
        raise ValueError("ka operator input support does not match kb output support")

    expected_kb = maxdiag(cv_diag, kb_provenance["stage_leakage_diagnostics"])
    observed_kb = validate_leakage_map(
        kb_meta.get("leakage_max"), require_all=True, context="kb cumulative leakage"
    )
    if observed_kb != expected_kb:
        raise ValueError("kb cumulative leakage does not equal max(cv, kb-stage)")

    expected_ka = maxdiag(observed_kb, ka_provenance["stage_leakage_diagnostics"])
    observed_ka = validate_leakage_map(
        ka_meta.get("leakage_max"), require_all=True, context="ka cumulative leakage"
    )
    if observed_ka != expected_ka:
        raise ValueError("ka cumulative leakage does not equal max(kb-cumulative, ka-stage)")

    return {
        "cv_output_support": int(cv_support),
        "kb_output_support": int(kb_support),
        "kb_operator_provenance": kb_provenance,
        "ka_operator_provenance": ka_provenance,
        "cv_to_kb_max_aggregation_verified": True,
        "kb_to_ka_max_aggregation_verified": True,
    }


def load_result(d,stage,ijk,freeze):
    stem=prefix(d,stage,ijk);path=stem.with_suffix(".npz");jpath=stem.with_suffix(".json")
    if not path.exists() or not jpath.exists():
        raise FileNotFoundError("Missing prior stage "+str(stem))
    meta=json.loads(jpath.read_text(encoding="utf-8"))
    if meta.get("freeze")!=freeze or meta.get("stage")!=stage or meta.get("channel")!=list(ijk):
        raise ValueError("Checkpoint metadata mismatch "+str(stem))
    if meta["sha256"]!=digest(path):raise ValueError("Checkpoint SHA256 mismatch "+str(stem))
    validate_leakage_map(meta.get("leakage_max"), require_all=False, context=str(stem))
    state=restore_state(path)
    if meta["support"]!=len(state) or not math.isclose(meta["norm"],norm(state),rel_tol=1e-12,abs_tol=1e-12):
        raise ValueError("Checkpoint support or norm mismatch "+str(stem))
    return state,meta
def save_result(d,stage,ijk,state,freeze,diag,extra=None):
    stem=prefix(d,stage,ijk);path=stem.with_suffix(".npz")
    persist_state(path,state)
    meta=dict(stage=stage,channel=list(ijk),freeze=freeze,
              sha256=digest(path),support=len(state),norm=norm(state),
              max_spin=max_spin(state),
              leakage_max=validate_leakage_map(diag, require_all=False, context=f"{stage}_{ijk}"),
              elapsed_s=round(time.perf_counter()-START,3))
    meta.update(extra or {})
    write_atomic(stem.with_suffix(".json"),lambda f:f.write((json.dumps(meta,indent=2,sort_keys=True)+"\n").encode()))
    print("CHECKPOINT",stage,ijk,"states",len(state),"norm",meta["norm"],"sha",meta["sha256"][:12],flush=True)
def maxdiag(a,b):
    left = validate_leakage_map(a, require_all=False, context="left stage diagnostics")
    right = validate_leakage_map(b, require_all=False, context="right stage diagnostics")
    merged = dict(left)
    for key, value in right.items():
        merged[key] = max(merged[key], value) if key in merged else value
    return merged
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage",choices=STAGES,required=True)
    p.add_argument("--source",type=int,default=1)
    p.add_argument("--column-index",type=int,default=0)
    p.add_argument("--a",type=int,default=0);p.add_argument("--b",type=int,default=2);p.add_argument("--c",type=int,default=3)
    p.add_argument("--coefficient",type=int,default=-1)
    p.add_argument("--channel",choices=["all"]+["".join(map(str,x)) for x in itertools.product(range(2),repeat=3)],default="all")
    p.add_argument("--checkpoint-dir",type=Path,required=True)
    p.add_argument("--json-output",type=Path)
    p.add_argument("--state-output",type=Path)
    args=p.parse_args()
    if args.source not in PW.VERT or len({args.a,args.b,args.c})!=3 or not all(v in PW.NEIG[args.source] for v in (args.a,args.b,args.c)):
        p.error("triple must be three distinct neighbors")
    if args.coefficient not in (-1,1):p.error("coefficient must be +/-1")
    basis=PW.basis_full_jhalf()
    if not 0<=args.column_index<len(basis):p.error("column-index out of 32")
    initial=basis[args.column_index]
    freeze=dict(source=args.source,column_index=args.column_index,input_key=repr(initial),
                input_Ks=list(initial[1]),triple=[args.a,args.b,args.c],epsilon=args.coefficient,
                JMAX2=JMAX2,TOL=TOL,checkpoint_schema=CHECKPOINT_SCHEMA,
                runtime=runtime_fingerprint(),modules=source_hashes())
    args.checkpoint_dir.mkdir(parents=True,exist_ok=True)
    LP.JMAX2=JMAX2
    restore,caches=LP.install_sine_cached_stack()
    try:
        if args.stage=="cv":
            psi=CV.gauss_to_covariant({initial:1+0j},args.source)
            for i,j,k in selected(args):
                out,leak=RAW.COMP.C_volume_component(psi,args.source,args.c,k,i,JMAX2)
                save_result(args.checkpoint_dir,"cv",(i,j,k),out,freeze,
                            {"CV_complete_basis_leakage":leak})
        elif args.stage in ("kb","ka"):
            prev="cv" if args.stage=="kb" else "kb"
            for i,j,k in selected(args):
                x,meta=load_result(args.checkpoint_dir,prev,(i,j,k),freeze)
                diag=meta["leakage_max"]
                exact_empty_input = not bool(x)
                if exact_empty_input:
                    out={}
                    # Exact empty input is a mathematical zero, not an absent
                    # measurement. Persist its provenance separately from the
                    # numerical value so the finish gate can audit the reason.
                    d={"outer_complete_basis_leakage":0.0,
                       "internal_volume_sector_leakage":0.0,
                       "complete_charge_basis_leakage":0.0}
                elif args.stage=="kb":
                    out,d=RAW.KCOMP.C_K_component(x,args.source,args.b,j,k,JMAX2)
                else:
                    out,d=RAW.KCOMP.C_K_component(x,args.source,args.a,i,j,JMAX2)
                d={("CK_"+field):float(value) for field,value in d.items()
                   if field in ("outer_complete_basis_leakage","internal_volume_sector_leakage","complete_charge_basis_leakage")}
                stage_diag = dict(d)
                diag=maxdiag(diag,stage_diag)
                before=len(out)
                if args.stage=="kb":
                    out={key:amp for key,amp in out.items() if key[2] in (0,2)}
                zero_fields = (["CK_outer_complete_basis_leakage",
                                "CK_internal_volume_sector_leakage",
                                "CK_complete_charge_basis_leakage"]
                               if exact_empty_input else [])
                save_result(args.checkpoint_dir,args.stage,(i,j,k),out,freeze,diag,
                            {"pruned_by_frozen_scalar_channel":before-len(out),
                             "operator_input_support":len(x),
                             "operator_stage_applied":not exact_empty_input,
                             "diagnostic_provenance_scope":"current_operator_stage",
                             "stage_leakage_diagnostics":stage_diag,
                             "exact_empty_input_zero_fields":zero_fields})
        else:
            total={};diag={};branch=[]
            for i,j,k in itertools.product(range(2),repeat=3):
                cv_state,cv_meta=load_result(args.checkpoint_dir,"cv",(i,j,k),freeze)
                kb_state,kb_meta=load_result(args.checkpoint_dir,"kb",(i,j,k),freeze)
                state,meta=load_result(args.checkpoint_dir,"ka",(i,j,k),freeze)
                chain_provenance = validate_checkpoint_chain(
                    cv_meta, kb_meta, meta, len(cv_state), len(kb_state)
                )
                LP.add(total,state) # exactly old pruning/TOL semantics
                diag=maxdiag(diag,meta["leakage_max"])
                branch.append({"ijk":[i,j,k],"gauss_projection_pending":True,
                              "support":len(state),"sha256":meta["sha256"],
                              **chain_provenance})
            diag = validate_leakage_map(diag, require_all=True, context="finish stage")
            if not all(np.isfinite([z.real, z.imag]).all() for z in total.values()):
                raise ValueError("finish stage covariant state contains non-finite amplitudes")
            gauss,accept2,reject2=LGA.project_scalar_gauss(total,args.source,TOL)
            if not all(np.isfinite([z.real, z.imag]).all() for z in gauss.values()):
                raise ValueError("finish stage Gauss state contains non-finite amplitudes")
            covnorm=norm(total);gnorm=norm(gauss)
            frac=1.0 if accept2+reject2<1e-30 else accept2/(accept2+reject2)
            physical=max(diag[k] for k in LEAK_FIELDS)
            checks=dict(complete_required_diagnostics=set(diag)==set(LEAK_FIELDS),
                  finite_covariant_norm=math.isfinite(covnorm),
                  finite_gauss_norm=math.isfinite(gnorm),
                  finite_covariant_amplitudes=all(np.isfinite([z.real, z.imag]).all() for z in total.values()),
                  finite_gauss_amplitudes=all(np.isfinite([z.real, z.imag]).all() for z in gauss.values()),
                  physical_basis_and_charge_leakage=physical<1e-8,
                  scalar_closure_fraction_or_exact_zero=frac>1-1e-10,
                  nonscalar_rejected_norm=math.sqrt(max(reject2,0))<1e-8,
                  single_HL_spin_wall=max_spin(gauss)<=JMAX2/2+1e-12)
            report=dict(status="staged genuine full-state ordered triple",
                  passed=all(checks.values()),freeze=freeze,source_node=args.source,
                  ordered_edges=[args.a,args.b,args.c],epsilon_coefficient=args.coefficient,
                  input_Ks=list(initial[1]),column_index=args.column_index,
                  covariant_support=len(total),covariant_norm=covnorm,
                  gauss_support=len(gauss),gauss_norm=gnorm,gauss_max_spin=max_spin(gauss),
                  scalar_closure_fraction=frac,nonscalar_rejected_norm=math.sqrt(max(reject2,0)),
                  physical_acceptance_max_leakage=physical,
                  complete_charge_basis_leakage=diag["CK_complete_charge_basis_leakage"],
                  leakage_diagnostics=diag,checks=checks,branches=branch,
                  is_full_Hermitian_Lorentzian=False,is_full_HDA=False)
            if not args.json_output or not args.state_output:p.error("finish requires json-output and state-output")
            args.state_output.parent.mkdir(parents=True,exist_ok=True)
            original=list(sorted(gauss.items(),key=lambda kv:repr(kv[0])))
            np.savez_compressed(args.state_output,
                spins=np.asarray([key[0] for key,v in original],np.int16).reshape(len(original),len(PW.EDGES)),
                Ks=np.asarray([key[1] for key,v in original],np.int16).reshape(len(original),len(PW.VERT)),
                amp=np.asarray([v for key,v in original],complex))
            args.json_output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
            print(json.dumps(report,indent=2),flush=True)
            if not report["passed"]:return 1
        print("stage complete",args.stage,"elapsed_s",round(time.perf_counter()-START,3),flush=True)
        return 0
    finally:
        restore()
START=time.perf_counter()
if __name__=="__main__":raise SystemExit(main())
