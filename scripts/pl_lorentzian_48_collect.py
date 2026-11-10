#!/usr/bin/env python3
"""Collect 24 forward + 24 adjoint PL Lorentzian terms and form physical S.

No science target enters this reduction.  Each worker already projects only the
exact scalar Gauss channel; this collector applies the preregistered PL epsilon
coefficients and forms S=-i(L-Ldagger)/2.
"""
from __future__ import annotations
import argparse,json,math,re,traceback
from pathlib import Path
import numpy as np
TOL=1e-10

def load_state(path):
    with np.load(path, allow_pickle=False) as d:
        spins=np.asarray(d['spins']);Ks=np.asarray(d['Ks']);amp=np.asarray(d['amp'])
    if spins.ndim!=2 or Ks.ndim!=2 or amp.ndim!=1:
        raise ValueError(f"malformed NPZ array ranks in {path}")
    if spins.shape[1]!=32 or Ks.shape[1]!=16 or not(len(spins)==len(Ks)==len(amp)):
        raise ValueError(f"malformed NPZ key dimensions or row alignment in {path}")
    if not np.issubdtype(spins.dtype,np.integer) or not np.issubdtype(Ks.dtype,np.integer):
        raise ValueError(f"spin and K arrays must be integer-valued in {path}")
    if not np.isfinite(amp.real).all() or not np.isfinite(amp.imag).all():
        raise ValueError(f"non-finite amplitudes in {path}")
    out={}
    for spin_row,k_row,value in zip(spins,Ks,amp):
        key=(tuple(map(int,spin_row)),tuple(map(int,k_row)))
        if key in out:raise ValueError(f"duplicate sparse-state key in {path}: {key}")
        out[key]=complex(value)
    return out
def add(dst,src,scale=1.0):
    for k,a in src.items():
        if not np.isfinite([a.real,a.imag]).all():raise ValueError(f"non-finite input amplitude for {k}")
        z=dst.get(k,0j)+scale*a
        if not np.isfinite([z.real,z.imag]).all():raise ValueError(f"non-finite accumulated amplitude for {k}")
        dst[k]=z
def prune(s,tol=TOL):
    if not all(np.isfinite([a.real,a.imag]).all() for a in s.values()):
        raise ValueError("refusing to prune a sparse state containing NaN or infinity")
    return {k:a for k,a in s.items() if abs(a)>tol}
def norm(s):return math.sqrt(sum(abs(a)**2 for a in s.values()))
def max_spin(s):return max((max(k[0]) for k in s),default=0)/2.0

def save_bundle(path,forward,adjoint,S,nedges=32,nverts=16):
    out={}
    for name,state in [('forward',forward),('adjoint',adjoint),('S',S)]:
        rows=sorted(state.items(),key=lambda kv:repr(kv[0]))
        if rows:
            out[name+'_spins']=np.asarray([k[0] for k,_ in rows],dtype=np.int16)
            out[name+'_Ks']=np.asarray([k[1] for k,_ in rows],dtype=np.int16)
            out[name+'_amp']=np.asarray([a for _,a in rows],dtype=np.complex128)
        else:
            out[name+'_spins']=np.zeros((0,nedges),dtype=np.int16);out[name+'_Ks']=np.zeros((0,nverts),dtype=np.int16);out[name+'_amp']=np.zeros((0,),dtype=np.complex128)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**out)

def run(root):
    metas=[]
    for p in Path(root).rglob('term_*.json'):
        d=json.loads(p.read_text())
        if d.get('passed') is not True:raise RuntimeError(f'worker must explicitly pass with boolean true: {p}')
        required_checks=('all_required_leakage_diagnostics_present_and_finite','physical_basis_and_charge_leakage')
        if not isinstance(d.get('checks'),dict) or any(d['checks'].get(k) is not True for k in required_checks):
            raise RuntimeError(f'missing/failed charge-complete diagnostics in {p}')
        for key in ('physical_acceptance_max_leakage','complete_charge_diagnostic','nonscalar_rejected_norm','scalar_closure_fraction'):
            value=d.get(key)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(float(value)) or float(value)<0.0:
                raise RuntimeError(f'invalid required numeric metadata {key} in {p}')
        if d['scalar_closure_fraction']>1.0+1e-10:
            raise RuntimeError(f'invalid scalar closure fraction in {p}')
        if d.get('exact_zero_ordered_term') not in (True,False):
            raise RuntimeError(f'missing exact-zero declaration in {p}')
        mode=d['mode'];idx=int(d['index'])
        if mode not in ('forward','adjoint') or not 0<=idx<24:
            raise RuntimeError(f'invalid worker identity in {p}')
        # The state paired with this metadata must have the exact same stem.
        # A previous implementation selected the first arbitrary *.npz in the
        # parent directory, which is unsafe in merged/materialized layouts.
        exact=p.with_suffix('.npz')
        if exact.exists():
            npz=exact
        else:
            cand=list(Path(root).rglob(f'term_{mode}_{idx}.npz'))
            if len(cand)!=1:raise RuntimeError(f'cannot uniquely resolve NPZ for {mode} {idx}: {cand[:8]}')
            npz=cand[0]
        metas.append((mode,idx,p,npz,d))
    if len(metas)!=48:raise RuntimeError(f'need 48 passed worker metadata files, got {len(metas)}')
    seen={(m,i) for m,i,_,_,_ in metas}
    expect={(m,i) for m in ('forward','adjoint') for i in range(24)}
    if seen!=expect:raise RuntimeError(f'worker orbit mismatch missing={sorted(expect-seen)} extra={sorted(seen-expect)}')

    sums={'forward':{},'adjoint':{}};rows=[]
    for mode,idx,p,npz,d in sorted(metas,key=lambda x:(x[0],x[1])):
        st=load_state(npz);coef=int(d['PL_epsilon_coefficient']);add(sums[mode],st,coef)
        rows.append({'mode':mode,'index':idx,'coefficient':coef,'ordered_target_nodes':d['ordered_target_nodes'],
                     'worker_support':len(st),'worker_norm':norm(st),'exact_zero_ordered_term':d['exact_zero_ordered_term'],
                     'max_leakage':d['physical_acceptance_max_leakage'],'scalar_closure_fraction':d['scalar_closure_fraction']})
    L=prune(sums['forward']);Ld=prune(sums['adjoint'])
    tmp={};add(tmp,L,+1);add(tmp,Ld,-1);S=prune({k:-0.5j*a for k,a in tmp.items()})

    seed=((1,)*32,(0,)*16);sdiag=S.get(seed,0j)
    seed_parity=sum(seed[0])%2
    wrong=sum(1 for k in S if sum(k[0])%2!=seed_parity)
    finite=all(np.isfinite([z.real,z.imag]).all() for st in (L,Ld,S) for z in st.values())
    max_leak=max(float(d['physical_acceptance_max_leakage']) for *_,d in metas)
    max_rej=max(float(d['nonscalar_rejected_norm']) for *_,d in metas)
    checks={
      'all_48_workers_loaded':True,'full_24x2_orbit_unique':seen==expect,
      'finite_combined_amplitudes':finite,'max_worker_leakage_below_1e-8':max_leak<1e-8,
      'max_worker_nonscalar_rejection_below_1e-8':max_rej<1e-8,
      'S_preserves_even_valence_seed_parity':wrong==0,
      'S_diagonal_matrix_element_real':abs(sdiag.imag)<1e-9,
      'single_L_engine_spin_wall':max(max_spin(L),max_spin(Ld),max_spin(S))<=3.5+1e-12,
    }
    return L,Ld,S,{
      'status':'exact 16-cell PL-S3 physical Hermitian Lorentzian node column',
      'passed':bool(all(checks.values())),'science_status':'AMPLITUDE_PRECURSOR_S_NODE0',
      'checks':checks,'source_node':0,'Jmax':3.5,'worker_count':48,
      'forward_exact_zero_terms':sum(r['exact_zero_ordered_term'] for r in rows if r['mode']=='forward'),
      'adjoint_exact_zero_terms':sum(r['exact_zero_ordered_term'] for r in rows if r['mode']=='adjoint'),
      'forward_support':len(L),'forward_norm':norm(L),'adjoint_support':len(Ld),'adjoint_norm':norm(Ld),
      'S_support':len(S),'S_norm':norm(S),'S_seed_diagonal_amplitude':[float(sdiag.real),float(sdiag.imag)],
      'S_is_nonzero_science_output':bool(norm(S)>1e-10),'wrong_seed_parity_outputs':wrong,
      'max_spin_forward':max_spin(L),'max_spin_adjoint':max_spin(Ld),'max_spin_S':max_spin(S),
      'max_worker_physical_leakage':max_leak,'max_worker_nonscalar_rejected_norm':max_rej,
      'terms':rows,
      'definition':'S=-i(L_raw-L_raw^dagger)/2 using the complete 24-term forward and adjoint PL epsilon worker-index orbits',
      'interpretation':'Physical Hermitian Lorentzian amplitude column on the independent 16-cell PL-S3 habitat. Nonzero magnitude is reported, not required for an infrastructure PASS.',
      'scope_note':'One source-node S column on the homogeneous seed. Collective W0 still requires the remaining source nodes, route action and target-independent depth-2 image/compression.'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--json-output',type=Path,required=True);p.add_argument('--state-output',type=Path,required=True);a=p.parse_args()
    try:L,Ld,S,out=run(a.root);code=0 if out['passed'] else 1
    except Exception as exc:
        L=Ld=S={};out={'status':'collector exception','passed':False,'science_status':'INFRASTRUCTURE_DIAGNOSTIC','error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
    a.json_output.parent.mkdir(parents=True,exist_ok=True);a.json_output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8');save_bundle(a.state_output,L,Ld,S);print(json.dumps(out,indent=2,sort_keys=True));return code
if __name__=='__main__':raise SystemExit(main())
