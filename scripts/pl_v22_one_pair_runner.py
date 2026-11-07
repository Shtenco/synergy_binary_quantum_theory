#!/usr/bin/env python3
"""Resumable ORIGINAL 32-path C1 runner; no synthetic physical amplitudes."""
import argparse,hashlib,json,subprocess,sys
from itertools import product
from pathlib import Path
from v21_stage_logic import operations
from pl_v22_stage_worker import load_charged

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def head(repo):return subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()

def accepted(prefix,stage,operation,commit,backend,previous,expected):
    p=Path(prefix);j=p.with_suffix('.json');n=p.with_suffix('.npz')
    if not j.exists() or not n.exists():return False
    try:
        m=json.loads(j.read_text())
        if not(m.get('passed') is True and m.get('stage')==stage
            and m.get('source_commit')==commit and m.get('volume_backend_sha256')==backend
            and m.get('guard',{}).get('stage_physical_basis_pass') is True):return False
        if tuple(m.get('operator_signature',()))!=tuple(operation):return False
        for k in ('source_node','omit','cycle','mode','Jmax2','input_state_sha256'):
            if m.get(k)!=expected[k]:return False
        if stage>1:
            if m.get('order')!=expected['order'] or m.get('aux')!=expected['aux']:return False
            prev=Path(previous)
            if m.get('predecessor_npz_sha256')!=sha(prev.with_suffix('.npz')):return False
            if m.get('predecessor_metadata_sha256')!=sha(prev.with_suffix('.json')):return False
        load_charged(n,m['npz_sha256'])
        return True
    except (OSError,ValueError,KeyError,TypeError):return False

def run(a):
    from pl_dual_complex import DualComplex,seed_16cell_boundary
    from pl_lorentzian_physical_pair_worker_v4 import pair_indices
    from pl_peter_weyl_euclidean import PLPeterWeylEuclidean
    repo=Path(__file__).resolve().parent.parent
    commit=head(repo)
    if a.commit and a.commit!=commit:raise ValueError('FROZEN_COMMIT_MISMATCH')
    D=DualComplex(seed_16cell_boundary());G=PLPeterWeylEuclidean(D)
    p,q,ip,iq=pair_indices(D,a.source,a.omit,a.cycle)
    targets={'p':tuple(D.neighbor[(a.source,r)] for r in p),
             'q':tuple(D.neighbor[(a.source,r)] for r in q)}
    base=a.out.resolve();base.mkdir(parents=True,exist_ok=True)
    backend=sha(repo/'scripts/tetrahedral_volume_backend.py')
    seed=((1,)*len(G.EDGES),(0,)*D.n_tets)
    seed_hash=hashlib.sha256(repr(seed).encode()).hexdigest()
    result={'source_commit':commit,'source_node':a.source,'omit':a.omit,'cycle':a.cycle,
       'Jmax2':7,'volume_backend_sha256':backend,'input_state_sha256':seed_hash,
       'paths':[],'independent_v2':{},'HDA_proved':False}
    count=0;first={}
    def command(cmd):
        nonlocal count
        subprocess.run(cmd,cwd=repo,check=True);count+=1
    def unfinished():
        result['new_operator_calls']=count
        result['completed_paths']=len(result['paths'])
        (base/'INCOMPLETE_MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
        return {'status':'INCOMPLETE_NOT_PHYSICS','passed':False,**result}
    def stage(n,path,spec,previous=None):
        expected={'source_node':a.source,'omit':a.omit,'cycle':a.cycle,'mode':spec['mode'],
                  'order':spec['order'],'aux':list(spec['aux']),'Jmax2':7,'input_state_sha256':seed_hash}
        if accepted(path,n,spec['ops'][n-1],commit,backend,previous,expected):return True
        if a.max_new and count>=a.max_new:return False
        cmd=[sys.executable,str(repo/'scripts/pl_v22_stage_worker.py'),
             '--stage',str(n),'--commit',commit,'--volume-backend',a.volume_backend,
             '--source-node',str(a.source),'--omit',str(a.omit),'--cycle',str(a.cycle),
             '--mode',spec['mode'],'--order',spec['order'],'--aux',*map(str,spec['aux']),
             '--output-prefix',str(path)]
        if previous:cmd+=['--previous-prefix',str(previous)]
        command(cmd)
        if not accepted(path,n,spec['ops'][n-1],commit,backend,previous,expected):
            raise ValueError('FAILED_RELOADING_ORIGINAL_CHECKPOINT')
        return True
    for mode,order in product(('forward','adjoint'),('p','q')):
        for aux in product((0,1),repeat=3):
            ops=operations(targets[order],mode,aux)
            spec={'mode':mode,'order':order,'aux':aux,'ops':ops}
            sig=tuple(ops[0])
            if sig not in first:first[sig]=base/'first'/('_'.join(map(str,sig)))
            one=first[sig]
            ident=f"{mode}_{order}_{''.join(map(str,aux))}"
            two=base/'paths'/ident/'second'
            three=base/'paths'/ident/'third'
            if not stage(1,one,spec):return unfinished()
            if not stage(2,two,spec,one):return unfinished()
            if not stage(3,three,spec,two):return unfinished()
            result['paths'].append({'stage1':str(one),'stage2':str(two),'stage3':str(three)})
    refs={'fp':('forward',ip),'fq':('forward',iq),
          'ap':('adjoint',ip),'aq':('adjoint',iq)}
    for label,(mode,index) in refs.items():
        prefix=base/'v2'/label; z=prefix.with_suffix('.npz');j=prefix.with_suffix('.json')
        valid=False
        if z.exists() and j.exists():
            try:
                m=json.loads(j.read_text())
                valid=(m.get('passed') is True and m.get('producer_kind')=='v22_independent_original_v2_worker'
                    and m.get('source_commit')==commit and m.get('volume_backend_sha256')==backend
                    and m.get('npz_sha256')==sha(z) and m.get('index')==index and m.get('mode')==mode)
            except (OSError,ValueError,KeyError):pass
        if not valid:
            if a.max_new and count>=a.max_new:return unfinished()
            command([sys.executable,str(repo/'scripts/pl_v22_independent_v2.py'),
                     '--mode',mode,'--index',str(index),'--source',str(a.source),
                     '--commit',commit,'--out',str(prefix)])
        result['independent_v2'][label]={'npz':str(z),'sha256':sha(z),
          'metadata_json':str(j),'metadata_sha256':sha(j)}
    manifest=base/'C1_MANIFEST.json'
    manifest.write_text(json.dumps(result,indent=2)+'\n')
    subprocess.run([sys.executable,str(repo/'scripts/pl_v22_pair_collector.py'),
       '--manifest',str(manifest),'--output-prefix',str(base/'C1_PAIR_RESULT')],cwd=repo,check=True)
    return {'status':'ONE_VERIFIED_V4_PAIR_NOT_HDA','passed':True,'new_calls':count,'paths':32}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=int,default=0)
    p.add_argument('--omit',type=int,choices=range(4),default=0)
    p.add_argument('--cycle',type=int,choices=range(3),default=0)
    p.add_argument('--volume-backend',default='original_tetrahedral')
    p.add_argument('--commit')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--max-new',type=int,default=0)
    a=p.parse_args()
    print(json.dumps(run(a),indent=2))
if __name__=='__main__':main()
