#!/usr/bin/env python3
"""One exact distributed PL-S3 Lorentzian ordered-triple worker.

For the canonical 16-cell dual complex compute either

  forward: Tr[C_a(K) C_b(K) C_c(V)] |Omega>

or the exact adjoint-ordered column implied by
C(K)_ij^dag=-C(K)_ji and C(V)_ij^dag=C(V)_ji,

  adjoint: sum_ijk C_c(V)_ik C_b(K)_kj C_a(K)_ji |Omega>.

The physical v1.2 Hermitian structural block is assembled only by the collector,
S=-i(L_raw-L_raw^dag)/2.  This worker never fits a coefficient and never drops a
representation using a GR target.  Exact-zero ordered terms are valid.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math,subprocess,sys,traceback
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_lorentzian_gauss_action_gate as LGA
from pl_dual_complex import DualComplex,seed_16cell_boundary
from pl_peter_weyl_euclidean import PLPeterWeylEuclidean
from pl_covariant_backend import install_pl_graph
from tetrahedral_volume_backend import install_tetrahedral_volume_backend

JMAX2=7
TOL=1e-10

def parity(base,perm):
    idx=[base.index(x) for x in perm]
    inv=sum(idx[i]>idx[j] for i in range(3) for j in range(i+1,3))
    return -1 if inv%2 else 1

def add(dst,src,scale=1.0,tol=1e-11):
    for k,a in src.items():
        z=dst.get(k,0j)+scale*a
        if abs(z)>tol:dst[k]=z
        elif k in dst:del dst[k]
def norm2(s):
    """Historical squared-norm API used by the V2 group worker."""
    return float(sum(abs(amplitude)**2 for amplitude in s.values()))

def state_norm(s):
    norm=0.0
    for amplitude in s.values():
        norm=math.hypot(norm,abs(amplitude))
    return norm

def projection_rejected_norm(state,source):
    norm=0.0
    for key,amplitude in state.items():
        if LGA.covariant_scalar_to_gauss_key(key,source) is None:
            norm=math.hypot(norm,abs(amplitude))
    return norm

def max_spin(s):return max((max(k[0]) for k in s),default=0)/2.0

def ordered_spec(D,source,index):
    if not 0<=index<24:raise ValueError('index must be 0..23')
    omit=index//6; pi=index%6
    base=tuple(r for r in range(4) if r!=omit)
    perm=tuple(itertools.permutations(base))[pi]
    targets=tuple(D.neighbor[(source,r)] for r in perm)
    coef=D.local_sign(source,omit)*parity(base,perm)
    return omit,base,perm,targets,int(coef)

def update_diag(dst,d):
    for name,val in d.items():
        if isinstance(val,(int,float)):dst[name]=max(dst.get(name,0.0),float(val))

def ordered_dagger_state(seed,source,a,b,c):
    """Exact component-adjoint ordering for one raw triple."""
    psi=LP.CV.gauss_to_covariant({seed:1+0j},source)
    total={};diag={
      'CV_complete_basis_leakage':0.0,
      'CK_outer_complete_basis_leakage':0.0,
      'CK_internal_volume_sector_leakage':0.0,
      'CK_complete_charge_basis_leakage':0.0,
    }
    for i,j,k in itertools.product(range(2),repeat=3):
        s1,d1=LP.RAW.KCOMP.C_K_component(psi,source,a,j,i,JMAX2)
        update_diag(diag,{
          'CK_outer_complete_basis_leakage':d1['outer_complete_basis_leakage'],
          'CK_internal_volume_sector_leakage':d1['internal_volume_sector_leakage'],
          'CK_complete_charge_basis_leakage':d1['complete_charge_basis_leakage']})
        if not s1:continue
        s2,d2=LP.RAW.KCOMP.C_K_component(s1,source,b,k,j,JMAX2)
        update_diag(diag,{
          'CK_outer_complete_basis_leakage':d2['outer_complete_basis_leakage'],
          'CK_internal_volume_sector_leakage':d2['internal_volume_sector_leakage'],
          'CK_complete_charge_basis_leakage':d2['complete_charge_basis_leakage']})
        # Final C(V) is rank 0+1, so only source J=0,1 can return to scalar J=0.
        s2={key:amp for key,amp in s2.items() if key[2] in (0,2)}
        if not s2:continue
        s3,leak=LP.RAW.COMP.C_volume_component(s2,source,c,i,k,JMAX2)
        diag['CV_complete_basis_leakage']=max(diag['CV_complete_basis_leakage'],float(leak))
        add(total,s3)
    return total,diag

def save_state(path,state,nedges,nverts):
    rows=sorted(state.items(),key=lambda kv:repr(kv[0]))
    if rows:
        spins=np.asarray([k[0] for k,_ in rows],dtype=np.int16)
        Ks=np.asarray([k[1] for k,_ in rows],dtype=np.int16)
        amp=np.asarray([a for _,a in rows],dtype=np.complex128)
    else:
        spins=np.zeros((0,nedges),dtype=np.int16);Ks=np.zeros((0,nverts),dtype=np.int16);amp=np.zeros((0,),dtype=np.complex128)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,spins=spins,Ks=Ks,amp=amp)

def run(index,mode='forward',source=0):
    D=DualComplex(seed_16cell_boundary());G=PLPeterWeylEuclidean(D)
    omit,base,perm,targets,coef=ordered_spec(D,source,index)
    seed=((1,)*len(G.EDGES),(0,)*D.n_tets)
    try:
        source_commit=subprocess.check_output(
            ['git','rev-parse','HEAD'],text=True,timeout=10).strip()
    except (OSError,subprocess.SubprocessError) as exc:
        raise ValueError('V2_SOURCE_COMMIT_NOT_VERIFIABLE') from exc
    if len(source_commit)!=40 or any(ch not in '0123456789abcdef' for ch in source_commit):
        raise ValueError('BAD_V2_SOURCE_COMMIT')
    dirty=subprocess.check_output(
        ['git','status','--porcelain','--untracked-files=no'],text=True,timeout=10)
    if dirty.strip():
        raise ValueError('V2_TRACKED_WORKTREE_NOT_CLEAN')
    backend_file=HERE/'tetrahedral_volume_backend.py'
    backend_sha=hashlib.sha256(backend_file.read_bytes()).hexdigest()
    oldJ=LP.JMAX2
    with install_tetrahedral_volume_backend():
      with install_pl_graph(G):
        LP.JMAX2=JMAX2
        restore,caches=LP.install_sine_cached_stack()
        try:
            if mode=='forward':cov,diag=LP.ordered_triple_state(seed,source,*targets)
            elif mode=='adjoint':cov,diag=ordered_dagger_state(seed,source,*targets)
            else:raise ValueError('mode must be forward or adjoint')
            gauss,accepted2,rejected2=LGA.project_scalar_gauss(cov,source,TOL)
            cache_info={name:{'hits':fn.cache_info().hits,'misses':fn.cache_info().misses,'currsize':fn.cache_info().currsize} for name,fn in caches.items()}
        finally:
            restore();LP.JMAX2=oldJ
    physical=max(float(diag.get('CV_complete_basis_leakage',0.0)),float(diag.get('CK_outer_complete_basis_leakage',0.0)),float(diag.get('CK_internal_volume_sector_leakage',0.0)))
    charge=float(diag.get('CK_complete_charge_basis_leakage',0.0))
    cn=state_norm(cov);gn=state_norm(gauss)
    rejected_norm=projection_rejected_norm(cov,source)
    scalar_fraction=1.0 if cn==0.0 else gn/cn
    checks={
      'finite_covariant_norm':math.isfinite(cn),'finite_gauss_norm':math.isfinite(gn),
      'finite_projection_diagnostics':all(math.isfinite(float(x)) for x in (accepted2,rejected2)),
      'finite_covariant_amplitudes':all(np.isfinite([z.real,z.imag]).all() for z in cov.values()),
      'finite_gauss_amplitudes':all(np.isfinite([z.real,z.imag]).all() for z in gauss.values()),
      'physical_basis_volume_leakage':physical<1e-8,
      'scalar_closure_fraction_or_exact_zero':scalar_fraction>1-1e-10,
      'nonscalar_rejected_norm':rejected_norm<1e-8,
      'single_L_spin_wall':max_spin(gauss)<=JMAX2/2+1e-12,
      'PL_orientation_coefficient':coef in (-1,1),
    }
    return gauss,{
      'status':'exact distributed PL-S3 Lorentzian ordered triple','passed':bool(all(checks.values())),
      'mode':mode,'index':index,'source_node':source,'omitted_local_slot':omit,
      'source_commit':source_commit,'volume_backend':'tetrahedral-volume-v1',
      'volume_backend_sha256':backend_sha,
      'base_local_slots':list(base),'permuted_local_slots':list(perm),'ordered_target_nodes':list(targets),
      'PL_epsilon_coefficient':coef,'Jmax':JMAX2/2,'input_key':repr(seed),
      'covariant_support':len(cov),'covariant_norm':cn,'gauss_support':len(gauss),'gauss_norm':gn,
      'gauss_max_spin':max_spin(gauss),
      'exact_zero_ordered_term':not any(z != 0j for z in cov.values()) and not any(z != 0j for z in gauss.values()),
      'scalar_closure_fraction':scalar_fraction,'nonscalar_rejected_norm':rejected_norm,
      'accepted2_legacy_diagnostic':float(accepted2),'rejected2_legacy_diagnostic':float(rejected2),
      'physical_acceptance_max_leakage':physical,'complete_charge_diagnostic':charge,
      'cache_info':cache_info,'checks':checks,'weighted_here':False,
      'orientation_note':'Coefficient = tetrahedral local orientation sign times permutation parity; on boundary-4-simplex this reduces to orientation[v] times the historical K5 epsilon orbit.',
      'scope_note':'One of 24 forward or 24 adjoint exact ordered terms. Physical S is formed only after both complete sums are collected.'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--index',type=int,required=True);p.add_argument('--mode',choices=('forward','adjoint'),required=True);p.add_argument('--source',type=int,default=0);p.add_argument('--json-output',type=Path,required=True);p.add_argument('--state-output',type=Path,required=True)
    a=p.parse_args()
    try:state,out=run(a.index,a.mode,a.source);code=0 if out['passed'] else 1
    except Exception as exc:
        state={};out={'status':'worker exception','passed':False,'mode':a.mode,'index':a.index,'error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
    D=DualComplex(seed_16cell_boundary());save_state(a.state_output,state,len(D.dual_edges()),D.n_tets)
    a.json_output.parent.mkdir(parents=True,exist_ok=True);a.json_output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps(out,indent=2,sort_keys=True));return code
if __name__=='__main__':raise SystemExit(main())
