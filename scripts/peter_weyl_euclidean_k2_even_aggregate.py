#!/usr/bin/env python3
"""Aggregate 25 reusable K2-even shards and rank-reveal the enlarged even carrier.

Builds the 800-column sparse map Q from q_(w,v,i)=H_w H_v b_i, its Gram
G2=Q^dag Q, the exact boundary-return block X=V0^dag Q, and the boundary-
orthogonal residual Gram

    G2_perp = G2 - X^dag X.

Because V0 is orthonormal, rank(G2_perp) is the number of genuinely new even
directions contributed beyond the 32-dimensional boundary.  Dense ambient
Peter-Weyl matrices are never built; only a sparse union-of-support matrix for
the actually generated 800 columns is used.
"""
from __future__ import annotations

import argparse, glob, hashlib, json, sys
from pathlib import Path
import numpy as np
from scipy import sparse

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW


def array_hash(A):
    A=np.ascontiguousarray(A,dtype=np.complex128)
    return hashlib.sha256(A.view(np.float64).tobytes()).hexdigest()


def key_of(row): return (tuple(int(x) for x in row['spins']),tuple(int(x) for x in row['K_labels']))

def run(paths,output_npz:Path|None=None):
    shards=[]
    for p in paths:
        d=json.loads(Path(p).read_text(encoding='utf-8'))
        if d.get('schema')!='BQG_EUCLIDEAN_K2_EVEN_SHARD_V1' or not d.get('passed',False): raise RuntimeError(f'invalid K2 shard {p}')
        shards.append((Path(p),d))
    if len(shards)!=25: raise RuntimeError(f'expected 25 K2 shards, got {len(shards)}')
    labels={(int(d['target_node']),int(d['source_node'])) for _,d in shards}
    if labels!={(w,v) for w in range(5) for v in range(5)}: raise RuntimeError('K2 shard coverage must be exactly target x source = 5x5')
    source_hashes={d.get('source_packet_sha256') for _,d in shards};domains={d.get('source_domain_label') for _,d in shards};jmax={float(d['Jmax']) for _,d in shards}
    if len(source_hashes)!=1 or len(domains)!=1 or jmax!={2.5}: raise RuntimeError('K2 shard provenance mismatch')
    shards.sort(key=lambda x:(int(x[1]['target_node']),int(x[1]['source_node'])))

    columns=[];column_labels=[];declared_boundary=[];allfinite=True;maxspin=0.0
    for _,d in shards:
        w=int(d['target_node']);v=int(d['source_node'])
        if len(d.get('states',[]))!=32: raise RuntimeError(f'shard {(w,v)} missing second-action states')
        by_i={int(s['input_index']):s for s in d['states']}
        if set(by_i)!=set(range(32)): raise RuntimeError(f'shard {(w,v)} input coverage mismatch')
        for i in range(32):
            s=by_i[i];col={}
            for r in s['state']:
                z=complex(float(r['amp'][0]),float(r['amp'][1]));allfinite &= bool(np.isfinite(z.real) and np.isfinite(z.imag));k=key_of(r);col[k]=col.get(k,0j)+z
            columns.append(col);column_labels.append((w,v,i));maxspin=max(maxspin,float(s['max_spin']))
            declared_boundary.append({int(r['boundary_index']):complex(float(r['amp'][0]),float(r['amp'][1])) for r in s.get('boundary_return',[])})
    if len(columns)!=800 or len(set(column_labels))!=800: raise RuntimeError('expected 800 unique labelled K2 columns')

    basis=PW.basis_full_jhalf()
    if len(basis)!=32: raise RuntimeError('frozen boundary basis unavailable')
    boundary_keys=list(basis)
    union=set(boundary_keys)
    for col in columns: union.update(col)
    keys=sorted(union);row_index={k:r for r,k in enumerate(keys)}
    rr=[];cc=[];vv=[]
    for c,col in enumerate(columns):
        for k,z in col.items(): rr.append(row_index[k]);cc.append(c);vv.append(z)
    A=sparse.csc_matrix((np.asarray(vv,dtype=np.complex128),(rr,cc)),shape=(len(keys),800))
    G=(A.getH()@A).toarray();G=.5*(G+G.conj().T)
    bidx=[row_index[k] for k in boundary_keys]
    X=A[bidx,:].toarray()

    # Cross-check sparse extraction against each shard's independently emitted boundary-return list.
    bx_err=0.0
    for c,decl in enumerate(declared_boundary):
        for bi in range(32): bx_err=max(bx_err,abs(X[bi,c]-decl.get(bi,0j)))

    Gperp=G-X.conj().T@X;Gperp=.5*(Gperp+Gperp.conj().T)
    ge=np.linalg.eigvalsh(G);pe=np.linalg.eigvalsh(Gperp)
    gscale=max(float(np.max(np.abs(ge))),1.0);pscale=max(float(np.max(np.abs(pe))),1.0)
    gtol=1e-10*gscale;ptol=1e-10*pscale
    rankG=int(np.sum(ge>gtol));rankP=int(np.sum(pe>ptol));nullG=800-rankG;nullP=800-rankP
    boundary_return_norms=np.sqrt(np.sum(np.abs(X)**2,axis=0));col_norms=np.sqrt(np.maximum(np.real(np.diag(G)),0.0))
    with np.errstate(divide='ignore',invalid='ignore'):
        frac=np.where(col_norms>0,boundary_return_norms/col_norms,0.0)

    hard={
        'exact_25_shards':len(shards)==25,
        'exact_800_labelled_columns':len(columns)==800 and len(set(column_labels))==800,
        'all_sparse_amplitudes_finite':bool(allfinite),
        'spin_cutoff_respected':maxspin<=2.5+1e-12,
        'G2_hermitian':float(np.linalg.norm(G-G.conj().T))<3e-9,
        'G2_positive_semidefinite':float(np.min(ge))>-3e-9*gscale,
        'G2_perp_hermitian':float(np.linalg.norm(Gperp-Gperp.conj().T))<3e-9,
        'G2_perp_positive_semidefinite':float(np.min(pe))>-3e-9*pscale,
        'boundary_return_reconstruction_matches_shards':bx_err<2e-12,
    }
    if output_npz is not None:
        output_npz=Path(output_npz);output_npz.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(output_npz,G2=G,X_boundary=X,G2_perp=Gperp,labels=np.asarray(column_labels,dtype=int))

    return {
        'schema':'BQG_EUCLIDEAN_K2_EVEN_SPAN_V1','passed':bool(all(hard.values())),
        'science_status':'MEASURED_EUCLIDEAN_K2_EVEN_CARRIER',
        'source_packet_sha256':next(iter(source_hashes)),'source_domain_label':next(iter(domains)),'Jmax':2.5,
        'boundary_dimension':32,'labelled_two_hit_columns':800,'ambient_generated_support_union_dimension':len(keys),
        'G2':{'dimension':800,'rank':rankG,'nullity':nullG,'rank_tolerance':gtol,'eigenvalue_min':float(ge[0]),'eigenvalue_max':float(ge[-1]),'trace':float(np.trace(G).real),'frobenius_norm':float(np.linalg.norm(G)),'hash':array_hash(G)},
        'boundary_return':{'matrix_shape':[32,800],'frobenius_norm':float(np.linalg.norm(X)),'max_column_fraction':float(np.max(frac)),'mean_column_fraction':float(np.mean(frac)),'reconstruction_max_abs_error':float(bx_err),'hash':array_hash(X)},
        'G2_perp':{'definition':'G2-X^dagger X','rank_new_even_directions':rankP,'nullity':nullP,'rank_tolerance':ptol,'eigenvalue_min':float(pe[0]),'eigenvalue_max':float(pe[-1]),'trace':float(np.trace(Gperp).real),'frobenius_norm':float(np.linalg.norm(Gperp)),'hash':array_hash(Gperp)},
        'total_even_carrier_dimension_through_two_hits':32+rankP,
        'max_spin_reached':maxspin,'hard_integrity_checks':hard,
        'history_relevance':{
            'first_enlarged_even_constraint_generated_carrier_measured':True,
            'identified_with_master_krylov_MV0':False,
            'mu2_computed':False,
            'physical_projector_emitted':False,
            'next_required_operation':'construct/evaluate the physical master on this even carrier (or obtain equivalent master moments) with the frozen Hermitian constraint convention and HDA certification',
        },
        'claim_boundary':'Rank-revealed two-forward-H Euclidean even carrier built from actual Peter-Weyl states. It is the first enlarged even constraint-generated habitat beyond V0, but it is not silently identified with M_E V0 and therefore is not mu2 or P_phys.'
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--shard',type=Path,action='append');ap.add_argument('--shard-glob');ap.add_argument('--output-json',type=Path,required=True);ap.add_argument('--output-npz',type=Path,required=True);a=ap.parse_args();paths=list(a.shard or [])
    if a.shard_glob: paths.extend(Path(x) for x in glob.glob(a.shard_glob,recursive=True))
    out=run(paths,a.output_npz);a.output_json.parent.mkdir(parents=True,exist_ok=True);a.output_json.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps(out,indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
