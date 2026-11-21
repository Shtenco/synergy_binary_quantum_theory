#!/usr/bin/env python3
"""Fail-closed compressed signed-BQG five-block Schur -> metric -> TT extractor.

Input is the small orthonormal effective-basis matrix produced upstream after
true P-subspace whitening and target-independent Q residual whitening.  The
first 30 orthonormal basis directions span the original labelled five-block
metric carrier; the remaining directions are Q.  The raw labelled P Gram is
kept separately so the Schur result can be returned to the frozen q coordinate
frame before applying the measured q->metric map.

No full microscopic Hilbert matrix is materialized here.  No constraint
spectral variable is renamed physical frequency.  The only possible output is
c_micro_spatial_BQG; c_BQG_IR remains null until the separate connected
physical projector/history -> Gamma_phys chain is completed.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
try:import s4_tt_six_wilson_predictor as WILSON
except Exception:WILSON=None

RTOL=1e-10
GEOM_TOL=2e-9
LEADING_TOL=2e-6
FIT_TOL=2e-8
HERM_TOL=5e-9
P_DIM=30


def json_default(x):
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,np.ndarray):return x.tolist()
    raise TypeError(f'Object of type {type(x).__name__} is not JSON serializable')


def frozen_metric_map():
    a=1/math.sqrt(12);b=1/math.sqrt(6)
    return np.asarray([[a,0,0,0,0,a],[0,a,0,0,a,0],[0,0,a,a,0,0],[0,0,b,-b,0,0],[0,b,0,0,-b,0],[b,0,0,0,0,-b]],float)

def hermitian_defect(A):return float(np.linalg.norm(A-A.conj().T)/max(np.linalg.norm(A),1e-300))
def blockdiag(A,n):return np.kron(np.eye(n),A)
def load_metadata(z):
    raw=z['metadata_json'];raw=raw.item() if np.ndim(raw)==0 else raw
    if isinstance(raw,bytes):raw=raw.decode('utf-8')
    return json.loads(str(raw))
def _pair(x):
    try:a,b=x;return [int(a),int(b)]
    except Exception:return None

def validate_provenance(meta):
    e=[]
    if meta.get('actual_bqg_operator') is not True:e.append('actual_bqg_operator must be true')
    if meta.get('synthetic') is not False:e.append('synthetic must be false')
    if meta.get('provenance_level')!='compressed_microscopic_constraint':e.append('provenance_level must be compressed_microscopic_constraint')
    if meta.get('physical_history_1pi') is not False:e.append('physical_history_1pi must be false')
    if meta.get('operator_family')!='frozen_signed_gravitational_constraint':e.append('wrong operator_family')
    if meta.get('raw_P_gram_preserved') is not True:e.append('raw_P_gram_preserved must be true')
    comps=set(meta.get('operator_components',[]))
    if not {'H_E_sine','S'}<=comps:e.append('operator components missing H_E_sine/S')
    if 'R_op' in comps or meta.get('route_operator_included') is not False:e.append('R_op forbidden')
    c=meta.get('operator_coefficients_exact',{})
    if _pair(c.get('H_E_sine'))!=[-2,3]:e.append('H_E coefficient must be -2/3')
    if _pair(c.get('S'))!=[-32,9]:e.append('S coefficient must be -32/9')
    if not str(meta.get('source_commit','')).strip():e.append('source_commit required')
    if not meta.get('regulator'):e.append('regulator required')
    if meta.get('basis_closure_complete') is not True:e.append('basis_closure_complete required')
    if meta.get('target_fitting_used') is not False:e.append('target_fitting_used must be false')
    return not e,e


def gram_sqrt(K):
    K=np.asarray(K,complex);K=.5*(K+K.conj().T);vals,U=np.linalg.eigh(K);scale=max(float(np.max(np.abs(vals))),1.0)
    if int(np.sum(vals>RTOL*scale))!=K.shape[0]:raise RuntimeError('raw P Gram rank deficient')
    root=(U*np.sqrt(vals))@U.conj().T
    return K,root,{'min_eigenvalue':float(vals.min()),'max_eigenvalue':float(vals.max()),'condition':float(vals.max()/vals.min())}


def schur_zero_energy_compressed(C,p_dim=P_DIM):
    C=np.asarray(C,complex)
    if C.ndim!=2 or C.shape[0]!=C.shape[1] or C.shape[0]<=p_dim:raise RuntimeError('C_basis must be square with non-empty Q')
    hd=hermitian_defect(C)
    if hd>HERM_TOL:raise RuntimeError(f'C_basis not Hermitian: {hd:.3e}')
    C=.5*(C+C.conj().T);A=C[:p_dim,:p_dim];B=C[:p_dim,p_dim:];D=C[p_dim:,p_dim:]
    vals,U=np.linalg.eigh(D);scale=max(float(np.max(np.abs(vals))),1.0);nz=np.abs(vals)>RTOL*scale;z=~nz
    coupled=0.0
    if np.any(z):
        coupled=float(np.linalg.norm(B@U[:,z])/max(np.linalg.norm(C),1e-300))
        if coupled>20*RTOL:raise RuntimeError(f'GAPLESS_COUPLED_Q_MODE_REQUIRES_PROMOTION: {int(np.sum(z))} Q zero modes, relative coupling={coupled:.3e}')
    inv=np.zeros_like(vals);inv[nz]=1/vals[nz];Dinv=(U*inv)@U.conj().T
    Ce=A-B@Dinv@B.conj().T;Ce=.5*(Ce+Ce.conj().T)
    av=np.abs(vals[nz]);gap=float(av.min()) if av.size else None;radius=float(av.max()) if av.size else None
    qres=D@(-Dinv@B.conj().T)+B.conj().T
    pres=A+B@(-Dinv@B.conj().T)-Ce
    return {'Ceff_orth':Ce,'q_dimension':C.shape[0]-p_dim,'q_eigenvalues':vals,'q_gap':gap,'q_spectral_radius':radius,
            'q_condition':None if gap in (None,0) else float(radius/gap),'zero_mode_count':int(np.sum(z)),'coupled_zero_mode_norm':coupled,
            'Schur_Q_residual':float(np.linalg.norm(qres)/max(np.linalg.norm(B),1e-300)),
            'Schur_P_residual':float(np.linalg.norm(pres)/max(np.linalg.norm(Ce),np.linalg.norm(A),1e-300))}


def canonicalize(Ce,K,pb,pc,positions,central):
    pb=np.asarray(pb,int).ravel();pc=np.asarray(pc,int).ravel();blocks=sorted(set(map(int,pb)))
    if len(blocks)!=5 or central not in blocks or pb.size!=30 or pc.size!=30:raise RuntimeError('invalid five-block P labels')
    if np.asarray(positions).shape!=(5,3):raise RuntimeError('positions must be 5x3 ordered by sorted block id')
    ob=[central]+[b for b in blocks if b!=central];order=[]
    for b in ob:
        for c in range(6):
            hit=np.flatnonzero((pb==b)&(pc==c))
            if hit.size!=1:raise RuntimeError(f'block {b} coord {c} not unique')
            order.append(int(hit[0]))
    return Ce[np.ix_(order,order)],K[np.ix_(order,order)],ob,blocks


def geometry(ob,blocks,positions):
    p={b:np.asarray(positions[i],float) for i,b in enumerate(blocks)};r0=p[ob[0]];offs=np.asarray([p[b]-r0 for b in ob[1:]],float);lens=np.linalg.norm(offs,axis=1)
    if np.min(lens)<=0:raise RuntimeError('zero neighbor displacement')
    a=float(np.mean(lens));u=offs/lens[:,None];sd=float(np.linalg.norm(u.sum(axis=0)));md=float(np.linalg.norm(u.T@u-(4/3)*np.eye(3)));ld=float(np.max(np.abs(lens/a-1)))
    return {'offsets':offs,'unit_normals':u,'a_star':a,'sum_defect':sd,'second_moment_defect':md,'equal_length_defect':ld,'passed':max(sd,md,ld)<GEOM_TOL}

def hvec(H):
    s=math.sqrt(2);return np.asarray([H[0,0],H[1,1],H[2,2],s*H[0,1],s*H[0,2],s*H[1,2]],float)
def tt_frame(n):
    if WILSON is None:raise RuntimeError('missing six-Wilson predictor')
    hp,hx=WILSON.tt_basis(n);return np.column_stack((hvec(hp),hvec(hx)))
def central_offsets(Kh):return Kh[:6,:6],[Kh[:6,6+6*i:12+6*i] for i in range(4)]
def taylor(K0,Knb,offs,n):
    n=np.asarray(n,float);n/=np.linalg.norm(n);mass=K0.copy();K2=np.zeros((6,6));K4=np.zeros((6,6));odd=0.0
    for A,r in zip(Knb,offs):
        A=np.asarray(A,float);S=A+A.T;D=A-A.T;x=float(n@r);mass+=S;K2+=-.5*S*x*x;K4+=(1/24)*S*x**4;odd=max(odd,float(np.linalg.norm(D)/max(np.linalg.norm(S),1e-300)))
    return mass,K2,K4,odd
DIRS={'100':np.array([1,0,0.]),'110':np.array([1,1,0.]),'111':np.array([1,1,1.]),'120':np.array([1,2,0.]),'generic':np.array([2,3,5.])}
EXTRACT_A=np.asarray([[1/6,0,0,0,0,0],[0,0,1/6,0,0,0],[5/96,1/48,0,1/96,1/24,1/96],[0,0,1/24,1/48,0,0],[1/81]*6,[341/3750,16/1875,0,17/1875,2/75,17/1875]],float)
def extract_six(rows):
    o=np.asarray([rows['100'][0,0],rows['100'][1,1],rows['110'][0,0],rows['110'][1,1],rows['111'][0,0],rows['120'][0,0]],float);return np.linalg.solve(EXTRACT_A,o)
def analyze_spatial(Kh,offs,a):
    K0,Knb=central_offsets(Kh);k2={};k4={};mass={};odd=0
    for name,n in DIRS.items():
        m,a2,a4,o=taylor(K0,Knb,offs,n);E=tt_frame(n);mass[name]=E.T@m@E;k2[name]=E.T@a2@E;k4[name]=E.T@a4@E;odd=max(odd,o)
    residues=np.asarray([np.trace(k2[k])/2 for k in DIRS]);Z=float(residues.mean())
    if abs(Z)<1e-14:raise RuntimeError('leading TT k2 residue is zero')
    ld=max(float(np.linalg.norm(k2[k]/Z-np.eye(2))) for k in DIRS);ms=abs(Z)/max(a*a,1e-300);md=max(float(np.linalg.norm(mass[k])/max(ms,1e-300)) for k in DIRS)
    dim={k:k4[k]/(Z*a*a) for k in DIRS};c=extract_six(dim);fd=0;pred={}
    for k,n in DIRS.items():
        p=np.asarray(WILSON.evaluate(c,n)['quartic_TT_matrix'],float);pred[k]=p;fd=max(fd,float(np.linalg.norm(dim[k]-p)/max(np.linalg.norm(dim[k]),1.0)))
    checks={'leading_k2_positive':Z>0,'leading_TT_isotropic':ld<LEADING_TOL,'TT_massless_at_k0':md<LEADING_TOL,'parity_even_reciprocal_transfer':odd<LEADING_TOL,'six_Wilson_fit_closes':fd<FIT_TOL}
    return {'Z2_spatial':Z,'leading_isotropy_defect':ld,'mass_defect':md,'reciprocity_odd_defect':odd,'wilson_fit_relative_defect':fd,'coefficients':c,'checks':checks,'passed':all(checks.values())}


def science_run(path,output=None):
    with np.load(path,allow_pickle=False) as z:
        req={'C_basis','P_gram','p_block','p_coord','block_positions','central_block','C00','metadata_json'};missing=sorted(req-set(z.files))
        if missing:raise RuntimeError(f'missing production arrays: {missing}')
        C=np.asarray(z['C_basis'],complex);K=np.asarray(z['P_gram'],complex);pb=np.asarray(z['p_block'],int);pc=np.asarray(z['p_coord'],int);positions=np.asarray(z['block_positions'],float);central=int(np.asarray(z['central_block']).item());C00=complex(np.asarray(z['C00']).item());meta=load_metadata(z);M=np.asarray(z['metric_map'],float) if 'metric_map' in z.files else frozen_metric_map()
    ok,errs=validate_provenance(meta)
    if not ok:raise RuntimeError('PRODUCTION_PROVENANCE_REJECTED: '+'; '.join(errs))
    if C.shape[0]<=30 or K.shape!=(30,30):raise RuntimeError('bad compressed dimensions')
    if M.shape!=(6,6) or np.linalg.matrix_rank(M,RTOL)!=6:raise RuntimeError('metric_map must be invertible 6x6')
    if abs(C00.imag)>1e-9:raise RuntimeError('C00 must be real')
    K,Kroot,kd=gram_sqrt(K);sch=schur_zero_energy_compressed(C,30)
    Ce_raw=Kroot@sch['Ceff_orth']@Kroot;Ce_raw=.5*(Ce_raw+Ce_raw.conj().T)
    Ce,Ke,ob,blocks=canonicalize(Ce_raw,K,pb,pc,positions,central)
    Kq=2*Ce.real-2*C00.real*Ke.real;T=blockdiag(np.linalg.inv(M),5);Kh=T.T@Kq@T
    geo=geometry(ob,blocks,positions)
    if not geo['passed']:raise RuntimeError('NEIGHBOUR_GEOMETRY_NOT_FROZEN_TETRAHEDRAL')
    spatial=analyze_spatial(Kh,geo['offsets'],geo['a_star']);micro=spatial['coefficients'].tolist() if spatial['passed'] else None
    out={'status':'BQG compressed signed-G five-block Schur-to-TT spatial extraction','science_status':'MICROSCOPIC_SIGNED_G_TT_SPATIAL_PRECURSOR_EXTRACTED' if spatial['passed'] else 'MICROSCOPIC_SIGNED_G_FAILS_IR_GUARDS','passed':bool(spatial['passed']),
         'source_metadata':meta,'carrier':{'dimension':30,**kd},'metric_map':M.tolist(),'metric_map_condition_number':float(np.linalg.cond(M)),
         'schur':{k:sch[k] for k in ('q_dimension','q_gap','q_spectral_radius','q_condition','zero_mode_count','coupled_zero_mode_norm','Schur_Q_residual','Schur_P_residual')},
         'geometry':{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in geo.items()},
         'spatial':{k:spatial[k] for k in ('Z2_spatial','leading_isotropy_defect','mass_defect','reciprocity_odd_defect','wilson_fit_relative_defect','checks')},
         'c_micro_spatial_BQG':micro,'c_BQG_IR':None,'physicalization_required_next':'actual constraints -> physical projector/history -> Z_phys[J_g] -> W_phys[J_g] -> Gamma_phys[g] -> Gamma_TT^(2)(omega,k)',
         'hard_scope_guard':'Constraint spatial precursor only. Physical c_BQG_IR remains null.'}
    if output is not None:output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(out,indent=2,default=json_default)+'\n',encoding='utf-8')
    return out


def selftest():
    if WILSON is None:raise RuntimeError('missing predictor')
    M=frozen_metric_map();normals=np.asarray([(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)],float)/math.sqrt(3)
    Kh=np.zeros((30,30))
    for b in range(5):Kh[6*b:6*b+6,6*b:6*b+6]=8*np.eye(6)
    for i in range(4):sl=slice(6+6*i,12+6*i);Kh[:6,sl]=-np.eye(6);Kh[sl,:6]=-np.eye(6)
    Bm=blockdiag(M,5);Kq_target=Bm.T@Kh@Bm
    L=np.eye(30);L[:-1,1:]+=0.08*np.eye(29);K=L.T@L;d=np.sqrt(np.diag(K));K=K/d[:,None]/d[None,:];K,Kroot,_=gram_sqrt(K);Kinv=np.linalg.inv(Kroot)
    Ce_raw=.5*Kq_target;Ce0=Kinv@Ce_raw@Kinv
    rng=np.random.default_rng(20260907);qdim=6;D=np.diag(np.linspace(2,5,qdim));X=rng.normal(size=(30,qdim))*.02;A=Ce0+X@np.linalg.inv(D)@X.T;C=np.block([[A,X],[X.T,D]])
    sch=schur_zero_energy_compressed(C,30);rec=Kroot@sch['Ceff_orth']@Kroot;err=float(np.linalg.norm(rec-Ce_raw)/np.linalg.norm(Ce_raw))
    Kq2=2*rec.real;T=blockdiag(np.linalg.inv(M),5);Kh2=T.T@Kq2@T;sp=analyze_spatial(Kh2,normals,1.0);exp=(-1/20)*WILSON.ISO+(1/18)*WILSON.Q4V;cerr=float(np.linalg.norm(sp['coefficients']-exp)/np.linalg.norm(exp))
    gap=False
    try:schur_zero_energy_compressed(np.asarray([[0.,.01],[.01,0.]],complex),1)
    except RuntimeError as e:gap='GAPLESS_COUPLED_Q_MODE_REQUIRES_PROMOTION' in str(e)
    meta={'actual_bqg_operator':True,'synthetic':False,'provenance_level':'compressed_microscopic_constraint','physical_history_1pi':False,'operator_family':'frozen_signed_gravitational_constraint','raw_P_gram_preserved':True,'operator_components':['H_E_sine','S'],'operator_coefficients_exact':{'H_E_sine':[-2,3],'S':[-32,9]},'route_operator_included':False,'source_commit':'SELFTEST','regulator':{'Jmax':'selftest'},'basis_closure_complete':True,'target_fitting_used':False}
    pok,pe=validate_provenance(meta);checks={'metric_map_cond_sqrt2':bool(abs(np.linalg.cond(M)-math.sqrt(2))<2e-12),'nonidentity_P_gram':bool(np.linalg.norm(K-np.eye(30))>1e-3),'compressed_Schur_recovers_raw_q':bool(err<5e-11),'gapless_coupled_Q_rejected':bool(gap),'provenance_accepts':bool(pok and not pe),'IR_guards_pass':bool(sp['passed']),'Wilson_known_answer':bool(cerr<2e-10)}
    return {'status':'compressed nonorthogonal-P known-answer test','science_status':'INFRASTRUCTURE_SELFTEST_NOT_BQG_EVIDENCE','passed':bool(all(checks.values())),'checks':checks,'schur_relative_error':err,'c_relative_error':cerr,'P_gram_condition':float(np.linalg.cond(K)),'c_BQG_IR':None,'hard_scope_guard':'Synthetic compressed-basis test only.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path);p.add_argument('--output',type=Path);p.add_argument('--selftest',action='store_true');a=p.parse_args()
    if a.selftest or a.input is None:o=selftest()
    else:
        try:o=science_run(a.input,a.output)
        except Exception as e:
            o={'status':'STOP','science_status':'MISSING_OR_INVALID_COMPRESSED_SIGNED_G_INPUT','passed':False,'c_micro_spatial_BQG':None,'c_BQG_IR':None,'error':str(e),'hard_scope_guard':'No Wilson coefficients emitted from incomplete/invalid compressed chain.'}
            if a.output is not None:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(o,indent=2,default=json_default)+'\n',encoding='utf-8')
    print(json.dumps(o,indent=2,default=json_default));return 0 if o.get('passed') else 2

if __name__=='__main__':raise SystemExit(main())