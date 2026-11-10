#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,pickle,tarfile,itertools
from pathlib import Path
import numpy as np
try:
    from .rank_matrix_guard import validate_full_column_matrix
except ImportError:
    from rank_matrix_guard import validate_full_column_matrix
import bqg_depth6_41 as Z
import bqg_depth6_generic_jucys_selector as J
import bqg_depth6_generic_master_row as M
import bqg_safe_boundary_master_gate as R
import bqg_s5_intertwiner_action_gate as S

TARGET={'32':(130903,2755),'311':(153455,2719),'221':(130503,2749),'2111':(103318,2712)}

def gauss_ok(sp):
    return all(R.allowed_k2_t(*R.local_spins(sp,v)) for v in R.VERT)

def geom(sp,v):
    sp=tuple(sp);Jm=max(sp)+1;o=set()
    for a,b in itertools.combinations(R.NEIG[v],2):
        es=[R.EIDX[tuple(sorted((v,a)))],R.EIDX[tuple(sorted((a,b)))],R.EIDX[tuple(sorted((b,v)))]]
        for ds in itertools.product((-1,1),repeat=3):
            z=list(sp);ok=True
            for e,dd in zip(es,ds):
                z[e]+=dd
                if z[e]<0 or z[e]>Jm:ok=False;break
            if ok:
                z=tuple(z)
                if gauss_ok(z):o.add(z)
    return o

def s2_bases_of_s5(rep):
    return tuple(sorted({M.canon2(S.permute_spins(rep,p))[0] for p in Z.S5}))

def support(rep,key):
    seed=J.CFG[key][0];outseed=M.flip(seed);out=set()
    for b in J.s3_bases_of_s5(tuple(rep)):
        for v in (0,1):
            for z in geom(b,v):
                q=Z.canon3(z)[0]
                if J.seed_basis(q,outseed).shape[1]:out.add((f'v{v}',q))
    for b in s2_bases_of_s5(tuple(rep)):
        for z in geom(b,2):
            q=M.canon2(z)[0]
            if M.basis2(q,outseed).shape[1]:out.add(('v2',q))
    return out

def extract(root,work,expected,prefix):
    root=Path(root);work=Path(work);work.mkdir(parents=True,exist_ok=True)
    ts=sorted(root.rglob(prefix+'*.tar.gz'))
    if len(ts)!=expected:raise SystemExit(f'expected {expected} tarballs got {len(ts)}')
    ds=[]
    for t in ts:
        d=work/t.stem.replace('.tar','');d.mkdir(exist_ok=True)
        with tarfile.open(t,'r:gz') as f:f.extractall(d)
        ds.append(d)
    return ds

def load_index(ds,key,expected):
    rec={};paths={};summ=[]
    for d in ds:
        ss=list(d.rglob('summary.json'))
        if len(ss)!=1:raise SystemExit(('summary',d,len(ss)))
        s=json.load(open(ss[0]));summ.append(s)
        if s['irrep']!=key or s['shell_states']!=264962 or s['s5_orbits']!=2757 or s['failed_blocks']!=0 or s['ok_blocks']!=s['assigned_blocks']:raise SystemExit(('bad summary',s))
        for r in s['records']:
            i=int(r['orbit_index'])
            if i in rec:raise SystemExit(('dup record',i))
            rec[i]=r
        for p in d.rglob('b*.pkl'):
            i=int(p.stem[1:])
            if i in paths:raise SystemExit(('dup map',i))
            paths[i]=p
    target,active=TARGET[key]
    if len(summ)!=expected or len(rec)!=active or len(paths)!=active or sum(int(r['m']) for r in rec.values())!=target:
        raise SystemExit(('coverage',len(summ),len(rec),len(paths),sum(int(r['m']) for r in rec.values()),target,active))
    return rec,paths

def certify(rec,paths,key,outfile):
    remaining=set(rec);sup={};occ={}
    for i,r in rec.items():
        ss=support(tuple(r['rep']),key);sup[i]=ss
        for q in ss:occ.setdefault(q,set()).add(i)
    rounds=[];cert={}
    while True:
        ready=[]
        for i in sorted(remaining):
            uq=[q for q in sup[i] if occ[q]=={i}]
            if not uq:continue
            with open(paths[i],'rb') as f:mp=pickle.load(f)
            mats=[mp[q] for q in uq if q in mp]
            del mp
            if not mats:continue
            A=np.vstack(mats);m=int(rec[i]['m']);A=validate_full_column_matrix(A,m,f"orbit {i}")
            if A.shape[0]<m:continue
            s=np.linalg.svd(A,compute_uv=False);smax=float(s[0]) if len(s) else 0.
            tol=max(A.shape)*np.finfo(float).eps*max(smax,1.)*100
            rank=int((s>tol).sum());sig=float(s[-1]) if len(s)>=m else 0.
            if rank==m and sig>max(tol,1e-10):ready.append((i,m,len(uq),A.shape[0],sig,tol))
        if not ready:break
        rounds.append({'round':len(rounds),'blocks':len(ready),'columns':sum(x[1] for x in ready)})
        for i,m,nq,nr,sig,tol in ready:cert[str(i)]={'m':m,'unique_rows':nq,'rows':nr,'sigma_min':sig,'rank_tol':tol}
        for i,*_ in ready:remaining.remove(i)
        for i,*_ in ready:
            for q in sup[i]:occ[q].discard(i)
        print('ROUND',rounds[-1],'remain',len(remaining),flush=True)
    target,active=TARGET[key]
    res={'status':'PASS' if not remaining else 'RESIDUAL_CORE','irrep':key,'target_blocks':active,'target_columns':target,'certified_blocks':len(cert),'certified_columns':sum(v['m'] for v in cert.values()),'remaining_blocks':len(remaining),'remaining_columns':sum(int(rec[i]['m']) for i in remaining),'rounds':rounds,'global_sigma_min':min((v['sigma_min'] for v in cert.values()),default=None),'remaining_ids':sorted(remaining),'proof_logic':'Exact one-row master identity B=A0+A1+3A2; conservative subgroup-geometric occupancy; each removed S5 orbit block has full-column-rank numeric master image on rows unique at removal.'}
    Path(outfile).write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps({k:res[k] for k in ['status','certified_blocks','certified_columns','remaining_blocks','remaining_columns','global_sigma_min']},indent=2),flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--irrep',choices=TARGET,required=True);ap.add_argument('--artifacts',required=True);ap.add_argument('--expected-shards',type=int,required=True);ap.add_argument('--prefix',required=True);ap.add_argument('--work',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    ds=extract(a.artifacts,a.work,a.expected_shards,a.prefix);rec,paths=load_index(ds,a.irrep,a.expected_shards);certify(rec,paths,a.irrep,a.out)
if __name__=='__main__':main()
