#!/usr/bin/env python3
"""Assemble frozen signed BQG gravity on the compressed five-block Krylov basis.

Production seam only; no new dynamics is introduced.

Input basis W=[Q_P,Q_Q] is orthonormal.  The first 30 vectors span the true
labelled five-block q carrier after P whitening; the remainder is the
orthogonal, target-independent Q complement.  The producer supplies

    E_basis = W^dag H_E^sine W
    S_basis = W^dag S W,  S=-i/2(L_raw-L_raw^dagger)

and this file forms only

    G_basis = -(2/3) E_basis -(32/9) S_basis.

The original 30x30 labelled-coordinate Gram P_gram is preserved for the
subsequent Schur/metric/TT extraction.  No dense microscopic Hilbert matrix is
required and no GR/TT/Wilson target is allowed in source selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

RTOL=1e-10
HERM_TOL=5e-9
GEOM_TOL=2e-9
GRAM_TOL=2e-9
LEAK_TOL=1e-9
E_COEFF=-2.0/3.0
S_COEFF=-32.0/9.0


def hermitian_defect(A):
    A=np.asarray(A,complex)
    return float(np.linalg.norm(A-A.conj().T)/max(np.linalg.norm(A),1e-300))


def sha256_file(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):
            h.update(chunk)
    return h.hexdigest()


def sha256_array(A)->str:
    A=np.ascontiguousarray(A)
    h=hashlib.sha256();h.update(str(A.dtype).encode());h.update(str(A.shape).encode());h.update(A.view(np.uint8))
    return h.hexdigest()


def load_metadata(z):
    raw=z['metadata_json'];raw=raw.item() if np.ndim(raw)==0 else raw
    if isinstance(raw,bytes):raw=raw.decode('utf-8')
    return json.loads(str(raw))


def validate_source_provenance(meta):
    errors=[]
    req={
        'synthetic':False,
        'target_fitting_used':False,
        'basis_closure_complete':True,
        'route_operator_included':False,
        'orthonormal_effective_basis':True,
        'p_first_in_effective_basis':True,
        'q_target_independent':True,
    }
    for k,v in req.items():
        if meta.get(k) is not v:errors.append(f'metadata.{k} must be {v!r}')
    if not str(meta.get('source_commit','')).strip():errors.append('metadata.source_commit is required')
    if not meta.get('regulator'):errors.append('metadata.regulator is required')
    if meta.get('basis_construction')!='P_whiten_then_Q_residual_whiten':errors.append('basis_construction must be P_whiten_then_Q_residual_whiten')
    comps=set(meta.get('operator_components',[]))
    if not {'H_E_sine','S'}<=comps:errors.append('operator_components must contain H_E_sine and S')
    if 'R_op' in comps:errors.append('R_op is forbidden in the signed gravitational TT precursor')
    defs=meta.get('component_definitions',{})
    if defs.get('H_E_sine') not in ('physical-sine Euclidean Hamiltonian','H_E_sine'):errors.append('component_definitions.H_E_sine is not frozen')
    if defs.get('S') not in ('-i/2*(L_raw-L_raw_dagger)','-i/2(L_raw-L_raw_dagger)'):errors.append('component_definitions.S must be Hermitian completion')
    leaks=meta.get('operator_compression_leakage_relative',{})
    for key in ('H_E_sine','S'):
        try:x=float(leaks[key])
        except Exception:
            errors.append(f'operator_compression_leakage_relative.{key} is required');continue
        if not math.isfinite(x) or x<0 or x>LEAK_TOL:errors.append(f'{key} compression leakage {x} exceeds {LEAK_TOL}')
    return not errors,errors


def validate_layout(pblock,pcoord,positions,central):
    pblock=np.asarray(pblock,int).ravel();pcoord=np.asarray(pcoord,int).ravel();positions=np.asarray(positions,float)
    if pblock.size!=30 or pcoord.size!=30:raise RuntimeError('five-block carrier requires 30 labelled q coordinates')
    blocks=sorted(set(map(int,pblock)))
    if len(blocks)!=5 or central not in blocks:raise RuntimeError('need exactly five blocks including central_block')
    for b in blocks:
        got=sorted(map(int,pcoord[pblock==b]))
        if got!=list(range(6)):raise RuntimeError(f'block {b} q labels invalid: {got}')
    if positions.shape!=(5,3):raise RuntimeError('block_positions must be 5x3 ordered by sorted block id')
    pos={b:positions[i] for i,b in enumerate(blocks)}
    offs=np.asarray([pos[b]-pos[central] for b in blocks if b!=central],float)
    lens=np.linalg.norm(offs,axis=1)
    if np.min(lens)<=0:raise RuntimeError('zero neighbor displacement')
    unit=offs/lens[:,None];a=float(np.mean(lens))
    d={
        'sum_defect':float(np.linalg.norm(np.sum(unit,axis=0))),
        'second_moment_defect':float(np.linalg.norm(unit.T@unit-(4/3)*np.eye(3))),
        'equal_length_defect':float(np.max(np.abs(lens/a-1))),
    }
    if max(d.values())>GEOM_TOL:raise RuntimeError('five-block shell is not tetrahedral: '+json.dumps(d))
    return {'blocks':blocks,'a_star':a,**d}


def validate_gram(K):
    K=np.asarray(K,complex)
    if K.shape!=(30,30):raise RuntimeError('P_gram must be 30x30')
    hd=hermitian_defect(K)
    if hd>HERM_TOL:raise RuntimeError(f'P_gram not Hermitian: {hd:.3e}')
    K=.5*(K+K.conj().T);vals=np.linalg.eigvalsh(K);scale=max(float(np.max(np.abs(vals))),1.0)
    rank=int(np.sum(vals>RTOL*scale))
    if rank!=30:raise RuntimeError(f'P_gram rank={rank}, expected 30')
    diagdef=float(np.max(np.abs(np.diag(K)-1)))
    if diagdef>GRAM_TOL:raise RuntimeError(f'P_gram diagonal is not unit-normalized q frame: {diagdef:.3e}')
    return K,{
        'rank':rank,'min_eigenvalue':float(vals.min()),'max_eigenvalue':float(vals.max()),
        'condition':float(vals.max()/vals.min()),'unit_diagonal_defect':diagdef,'Hermiticity_defect':hd,
    }


def assemble(input_path:Path,output_path:Path,summary_path:Path|None=None):
    with np.load(input_path,allow_pickle=False) as z:
        required={'E_basis','S_basis','P_gram','p_block','p_coord','block_positions','central_block','C00_E','C00_S','metadata_json'}
        missing=sorted(required-set(z.files))
        if missing:raise RuntimeError(f'missing source arrays: {missing}')
        E=np.asarray(z['E_basis'],complex);S=np.asarray(z['S_basis'],complex);Kraw=np.asarray(z['P_gram'],complex)
        pb=np.asarray(z['p_block'],int);pc=np.asarray(z['p_coord'],int);positions=np.asarray(z['block_positions'],float)
        central=int(np.asarray(z['central_block']).item());c00e=complex(np.asarray(z['C00_E']).item());c00s=complex(np.asarray(z['C00_S']).item())
        meta=load_metadata(z);metric_map=np.asarray(z['metric_map'],float) if 'metric_map' in z.files else None
        basis_ids=np.asarray(z['basis_ids']).astype(str) if 'basis_ids' in z.files else None
    ok,errs=validate_source_provenance(meta)
    if not ok:raise RuntimeError('SOURCE_PROVENANCE_REJECTED: '+'; '.join(errs))
    if E.ndim!=2 or E.shape[0]!=E.shape[1] or S.shape!=E.shape:raise RuntimeError('E_basis/S_basis must be same square shape')
    n=E.shape[0]
    if n<=30:raise RuntimeError('effective basis must include non-empty Q complement beyond 30 P directions')
    if basis_ids is not None and (basis_ids.shape!=(n,) or len(set(map(str,basis_ids)))!=n):raise RuntimeError('basis_ids must be N unique ids')
    ed=hermitian_defect(E);sd=hermitian_defect(S)
    if ed>HERM_TOL:raise RuntimeError(f'E_basis not Hermitian: {ed:.3e}')
    if sd>HERM_TOL:raise RuntimeError(f'S_basis not Hermitian: {sd:.3e}')
    if abs(c00e.imag)>HERM_TOL or abs(c00s.imag)>HERM_TOL:raise RuntimeError('C00_E/C00_S must be real')
    K,gdiag=validate_gram(Kraw);layout=validate_layout(pb,pc,positions,central)
    if metric_map is not None and (metric_map.shape!=(6,6) or np.linalg.matrix_rank(metric_map,RTOL)!=6):raise RuntimeError('metric_map must be invertible 6x6')
    G=E_COEFF*E+S_COEFF*S;G=.5*(G+G.conj().T);gd=hermitian_defect(G)
    C00=E_COEFF*c00e.real+S_COEFF*c00s.real
    outmeta={
        'actual_bqg_operator':True,'synthetic':False,'provenance_level':'compressed_microscopic_constraint','physical_history_1pi':False,
        'operator_family':'frozen_signed_gravitational_constraint','operator_components':['H_E_sine','S'],
        'operator_coefficients_exact':{'H_E_sine':[-2,3],'S':[-32,9]},'operator_formula':'G=-(2/3)H_E_sine-(32/9)S',
        'route_operator_included':False,'source_commit':meta['source_commit'],'regulator':meta['regulator'],'basis_closure_complete':True,
        'target_fitting_used':False,'effective_basis_order':'first 30 = whitened span(raw labelled P); remainder = orthogonal target-independent Q',
        'raw_P_gram_preserved':True,'source_bundle_sha256':sha256_file(input_path),'E_basis_sha256':sha256_array(E),'S_basis_sha256':sha256_array(S),
        'G_basis_sha256':sha256_array(G),'P_gram_sha256':sha256_array(K),'source_component_metadata':meta,
    }
    kw=dict(C_basis=G,P_gram=K,p_block=pb,p_coord=pc,block_positions=positions,central_block=np.asarray(central,int),C00=np.asarray(C00,float),metadata_json=np.asarray(json.dumps(outmeta,sort_keys=True)))
    if metric_map is not None:kw['metric_map']=metric_map
    if basis_ids is not None:kw['basis_ids']=basis_ids
    output_path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(output_path,**kw)
    out={
        'status':'assembled frozen signed G on compressed five-block Krylov basis','science_status':'COMPRESSED_SIGNED_G_READY_FOR_SCHUR','passed':True,
        'effective_dimension':n,'P_dimension':30,'Q_dimension':n-30,'E_coefficient':E_COEFF,'S_coefficient':S_COEFF,
        'E_Hermiticity_defect':ed,'S_Hermiticity_defect':sd,'G_Hermiticity_defect':gd,'C00_signed_G':float(C00),
        'P_gram':gdiag,'geometry':layout,'source_bundle_sha256':outmeta['source_bundle_sha256'],'G_basis_sha256':outmeta['G_basis_sha256'],
        'output':str(output_path),'hard_scope_guard':'Compressed microscopic signed constraint only; not Z_phys, Gamma_phys, physical omega, or c_BQG_IR.',
    }
    if summary_path is not None:
        summary_path.parent.mkdir(parents=True,exist_ok=True);summary_path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    return out


def selftest():
    rng=np.random.default_rng(20260907);n=36
    X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));E=.5*(X+X.conj().T)
    Y=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));S=.5*(Y+Y.conj().T)
    L=np.eye(30);L[:-1,1:]+=0.015*np.eye(29);K=L.T@L;d=np.sqrt(np.diag(K));K=K/d[:,None]/d[None,:]
    pb=np.repeat(np.arange(5),6);pc=np.tile(np.arange(6),5)
    normals=np.asarray([(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)],float)/math.sqrt(3);pos=np.vstack([np.zeros(3),normals])
    meta={'synthetic':False,'target_fitting_used':False,'basis_closure_complete':True,'source_commit':'SELFTEST','regulator':{'Jmax':'selftest'},
          'operator_components':['H_E_sine','S'],'route_operator_included':False,'component_definitions':{'H_E_sine':'H_E_sine','S':'-i/2(L_raw-L_raw_dagger)'},
          'orthonormal_effective_basis':True,'p_first_in_effective_basis':True,'q_target_independent':True,'basis_construction':'P_whiten_then_Q_residual_whiten',
          'operator_compression_leakage_relative':{'H_E_sine':0.0,'S':0.0}}
    ok,errs=validate_source_provenance(meta);Kg,gd=validate_gram(K);layout=validate_layout(pb,pc,pos,0);G=E_COEFF*E+S_COEFF*S
    bad=dict(meta);bad['route_operator_included']=True;bad['operator_components']=['H_E_sine','S','R_op'];bok,_=validate_source_provenance(bad)
    checks={
        'valid_provenance':bool(ok and not errs),
        'route_rejected':bool(not bok),
        'exact_signed_coefficients':bool(np.linalg.norm(G-(-(2/3)*E-(32/9)*S))<1e-12*np.linalg.norm(G)),
        'nonidentity_raw_P_gram':bool(np.linalg.norm(Kg-np.eye(30))>1e-5),
        'raw_P_rank30':bool(gd['rank']==30),
        'tetrahedral_layout':bool(max(layout[k] for k in ('sum_defect','second_moment_defect','equal_length_defect'))<GEOM_TOL),
        'Hermitian_components':bool(hermitian_defect(E)<HERM_TOL and hermitian_defect(S)<HERM_TOL),
    }
    return {'status':'compressed signed-G assembler selftest','science_status':'INFRASTRUCTURE_SELFTEST_NOT_BQG_EVIDENCE','passed':bool(all(checks.values())),
            'checks':checks,'P_gram_condition':float(gd['condition']),'hard_scope_guard':'Synthetic compressed-basis algebra test only.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path);p.add_argument('--output',type=Path);p.add_argument('--summary',type=Path);p.add_argument('--selftest',action='store_true');a=p.parse_args()
    if a.selftest:o=selftest()
    else:
        if a.input is None or a.output is None:p.error('--input and --output required unless --selftest')
        try:o=assemble(a.input,a.output,a.summary)
        except Exception as e:
            o={'status':'STOP','science_status':'MISSING_OR_INVALID_COMPRESSED_OPERATOR_SOURCE','passed':False,'error':str(e),'hard_scope_guard':'No signed matrix emitted from incomplete/fitted/route-mixed/leaky input.'}
            if a.summary is not None:
                a.summary.parent.mkdir(parents=True,exist_ok=True);a.summary.write_text(json.dumps(o,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(o,indent=2));return 0 if o.get('passed') else 2


if __name__=='__main__':raise SystemExit(main())
