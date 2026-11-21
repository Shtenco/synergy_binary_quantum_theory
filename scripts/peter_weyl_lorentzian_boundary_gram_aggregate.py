#!/usr/bin/env python3
"""Aggregate 32 complete raw Lorentzian columns at one K5 source node.

Inputs are BQG_LORENTZIAN_FULL_COLUMN_DAG_V1 artifacts.  The gate constructs
complete outgoing-state Gram M_L, the direct logical-return matrix R, audits the
source-stabilizer S4 commutant, and optionally checks the preregistered parity
orthogonality X_EL=0 against a reusable Euclidean boundary packet.

Symmetry is measured, never assumed.  A failed S4 commutant test does not erase
the measured Gram, but it explicitly sets symmetry_reconstruction_allowed=false.
"""
from __future__ import annotations

import argparse, glob, json, math, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
import bqg_constraint_master_assembler_gate as MASTER
import k5_logical_source_stabilizer_s4_gate as S4
import logical_s4_twirl_gate as LS4
import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL

TOL=2e-9


def decode_covariant(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']),tuple(int(x) for x in r['Kother']),int(r['J2']),int(r['M2']),int(r['K12']),int(r['K34']))
        out[key]=out.get(key,0j)+complex(float(r['amp'][0]),float(r['amp'][1]))
    return out


def decode_gauss(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']),tuple(int(x) for x in r['K_labels']))
        out[key]=out.get(key,0j)+complex(float(r['amp'][0]),float(r['amp'][1]))
    return out


def sparse_norm(s): return math.sqrt(sum(abs(z)**2 for z in s.values()))


def cross(A,B):
    n=len(A);m=len(B);X=np.zeros((n,m),complex)
    for i,a in enumerate(A):
        for j,b in enumerate(B):
            if len(a)>len(b): X[i,j]=np.conj(sum(np.conj(z)*a.get(k,0j) for k,z in b.items()))
            else: X[i,j]=sum(np.conj(z)*b.get(k,0j) for k,z in a.items())
    return X


def perm_sign_on_neighbors(g):
    p=[g[x] for x in (1,2,3,4)]
    inv=sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
    return -1 if inv%2 else 1


def run(column_paths,e_packet_dir=None):
    raw=[]
    for p in column_paths:
        d=json.loads(Path(p).read_text(encoding='utf-8'))
        if d.get('schema')!='BQG_LORENTZIAN_FULL_COLUMN_DAG_V1' or not d.get('passed',False): raise RuntimeError(f'invalid L column {p}')
        raw.append((Path(p),d))
    if len(raw)!=32: raise RuntimeError(f'expected 32 L columns, got {len(raw)}')
    idxs=[int(d['input_logical_basis_index']) for _,d in raw]
    if set(idxs)!=set(range(32)) or len(set(idxs))!=32: raise RuntimeError('L input coverage must be exactly 0..31')
    raw.sort(key=lambda x:int(x[1]['input_logical_basis_index']))
    sources={int(d['source_node']) for _,d in raw}; jmax={float(d['Jmax']) for _,d in raw}; hab={d['habitat_hash'] for _,d in raw}; dom={d['boundary_domain_hash'] for _,d in raw}; conv={d['convention_hash'] for _,d in raw}
    if len(sources)!=1 or len(jmax)!=1 or len(hab)!=1 or len(dom)!=1 or len(conv)!=1: raise RuntimeError('L column provenance mismatch')
    source=next(iter(sources))

    Lcov=[];Lgauss=[];isometry=[];R=np.zeros((32,32),complex)
    for col,(_,d) in enumerate(raw):
        s=decode_covariant(d['state']);g,mapdiag=FULL.project_covariant_J0_to_gauss(s,source)
        if mapdiag['invalid_J0_covariant_keys'] or mapdiag['mapping_collisions']: raise RuntimeError(f'bad J0 reverse map in column {col}')
        Lcov.append(s);Lgauss.append(g)
        nc=sparse_norm(s);ng=sparse_norm(g);isometry.append(abs(nc-ng)/max(nc,ng,1e-300))
        for row in d['logical_return']['nonzero_amplitudes']:
            R[int(row['logical_basis_index']),col]=complex(float(row['amp'][0]),float(row['amp'][1]))

    ML_cov=MASTER.gram(Lcov);ML_g=MASTER.gram(Lgauss)
    gram_map_rel=float(np.linalg.norm(ML_cov-ML_g)/max(np.linalg.norm(ML_cov),1e-300))
    ML=.5*(ML_cov+ML_cov.conj().T); audit=MASTER.spectral_audit(ML);ev=np.asarray(audit['eigenvalues'],float);scale=max(float(np.max(np.abs(ev))),1.0)

    # Raw logical-return diagnostics frozen before choosing a physical Hermitian convention.
    RH=.5*(R+R.conj().T);RA=(R-R.conj().T)/(2j);rn=max(float(np.linalg.norm(R)),1e-300)
    logical_diag={
        'frobenius_norm':float(np.linalg.norm(R)),
        'hermitian_defect_relative':float(np.linalg.norm(R-R.conj().T)/rn),
        'antihermitian_defect_relative':float(np.linalg.norm(R+R.conj().T)/rn),
        'R_H_frobenius_norm':float(np.linalg.norm(RH)),
        'R_A_frobenius_norm':float(np.linalg.norm(RA)),
    }

    # Source-stabilizer S4: Gram must commute whether raw H_L transforms as
    # scalar or sign-character. Direct R distinguishes the two hypotheses.
    if source!=0:
        symmetry={'evaluated':False,'reason':'current exact source-stabilizer implementation is frozen for source 0','reconstruction_allowed':False}
    else:
        lb=LS4.singlet_basis();perms=[(0,)+tuple(p) for p in __import__('itertools').permutations((1,2,3,4))]
        reps={g:S4.global_U(g,lb) for g in perms}; mln=max(float(np.linalg.norm(ML)),1e-300)
        gram_err=max(float(np.linalg.norm(U@ML@U.conj().T-ML)/mln) for U in reps.values())
        ordinary=max(float(np.linalg.norm(U@R@U.conj().T-R)/rn) for U in reps.values()) if np.linalg.norm(R)>0 else 0.0
        signerr=max(float(np.linalg.norm(U@R@U.conj().T-perm_sign_on_neighbors(g)*R)/rn) for g,U in reps.items()) if np.linalg.norm(R)>0 else 0.0
        symmetry={'evaluated':True,'gram_commutant_max_relative_error':gram_err,'raw_R_scalar_character_max_relative_error':ordinary,'raw_R_sign_character_max_relative_error':signerr,'gram_commutant_certified':gram_err<2e-8,'reconstruction_allowed':gram_err<2e-8,'note':'Gram commutant is necessary for symmetry reconstruction; held-out cross-source direct columns remain mandatory before replacing missing nodes.'}

    parity=None
    if e_packet_dir is not None:
        e_packet_dir=Path(e_packet_dir);E=[]
        for i in range(32):
            p=e_packet_dir/'columns'/f'E_node{source}_input{i:02d}.json'
            d=json.loads(p.read_text(encoding='utf-8'));E.append(decode_gauss(d['complete_gauss_outgoing_column']['state']))
        X=cross(E,Lgauss);xnorm=float(np.linalg.norm(X));den=max(math.sqrt(float(np.linalg.norm(MASTER.gram(E)))*float(np.linalg.norm(ML_g))),1e-300)
        parity={'X_EL_frobenius_norm':xnorm,'normalized_mixed_block':xnorm/den,'certified_zero':xnorm/den<2e-10}

    integrity={
        'exact_32_input_coverage':idxs==list(range(32)),
        'common_source_and_provenance':len(sources)==len(jmax)==len(hab)==len(dom)==len(conv)==1,
        'covariant_to_gauss_norm_isometry':max(isometry)<2e-10,
        'covariant_and_gauss_grams_agree':gram_map_rel<2e-10,
        'M_L_hermitian':float(np.linalg.norm(ML_cov-ML_cov.conj().T))<TOL,
        'M_L_positive_semidefinite':float(np.min(ev))>-TOL*scale,
    }
    if parity is not None: integrity['X_EL_zero_by_parity']=bool(parity['certified_zero'])
    return {
        'schema':'BQG_LORENTZIAN_BOUNDARY_GRAM_V1','passed':bool(all(integrity.values())),
        'source_node':source,'Jmax':next(iter(jmax)),'habitat_hash':next(iter(hab)),'boundary_domain_hash':next(iter(dom)),'convention_hash':next(iter(conv)),
        'M_L':{'rank':audit['rank'],'nullity':audit['nullity'],'rank_tolerance':audit['rank_tolerance'],'eigenvalue_min':float(np.min(ev)),'eigenvalue_max':float(np.max(ev)),'smallest_positive':audit['smallest_positive'],'condition_number_on_support':audit['condition_number_on_support'],'trace':float(np.trace(ML).real),'frobenius_norm':float(np.linalg.norm(ML)),'hash':MASTER.hash_arrays(ML)},
        'gauss_map':{'max_relative_norm_error':max(isometry),'gram_relative_error':gram_map_rel},
        'logical_return_matrix':logical_diag,'source_stabilizer_S4':symmetry,'parity_mixed_block':parity,'integrity_checks':integrity,
        'symmetry_reconstruction_allowed':bool(symmetry.get('reconstruction_allowed',False)),
        'claim_boundary':'Measured one-node raw Lorentzian boundary Gram and logical-return diagnostics. Even if symmetry passes, held-out direct cross-source columns are required before reconstructing missing node data. This is not the enlarged physical projector or HDA certificate.'
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--column',type=Path,action='append');ap.add_argument('--column-glob');ap.add_argument('--euclidean-packet-dir',type=Path);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    paths=list(a.column or [])
    if a.column_glob: paths.extend(Path(x) for x in glob.glob(a.column_glob))
    out=run(paths,a.euclidean_packet_dir);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps(out,indent=2));return 0 if out['passed'] else 1

if __name__=='__main__': raise SystemExit(main())
