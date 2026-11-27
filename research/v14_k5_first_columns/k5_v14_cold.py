#!/usr/bin/env python3
"""V14: cold K5 direct columns and selected independent reverse matrix elements.
Always run against canonical scripts/k5_peter_weyl_safe_hda_column.py.
Not a physical projector or a continuum calculation.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
SCRIPTS=HERE.parent.parent/'scripts'
sys.path.insert(0,str(SCRIPTS))
import k5_peter_weyl_safe_hda_column as K5

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def encode(k,a):return {'spins':list(k[0]),'Ks':list(k[1]),'re':float(a.real).hex(),'im':float(a.imag).hex()}
def decode(e):return ((tuple(e['spins']),tuple(e['Ks'])),complex(float.fromhex(e['re']),float.fromhex(e['im'])))
def atomic_json(p,o):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    t=p.with_suffix('.tmp');t.write_text(json.dumps(o,indent=2)+'\n');t.replace(p)

def run(folder:Path,jmax2:int,source_index:int,vertices:list[int],reverse_per_vertex:int):
    folder.mkdir(parents=True,exist_ok=True)
    source=K5.basis_full_jhalf()[source_index]
    meta={'sha256_upstream':sha(K5.__file__),'jmax2':jmax2,
          'source_index':source_index,'source_spins':list(source[0]),'source_Ks':list(source[1])}
    seen=set();weight=0.;by_vertex=[];reverse=[];changed=0.;nrows=0
    for v in vertices:
        f=folder/f'direct_v{v}_s{source_index}_wall{jmax2}.json'
        if f.exists():
            saved=json.loads(f.read_text())
            if saved['meta']!=meta:raise RuntimeError('stale provenance: '+str(f))
            out=dict(decode(a) for a in saved['rows'])
        else:
            t0=time.perf_counter()
            out=K5.prune_state(K5.apply_H_cached_state({source:1+0j},v,jmax2),1e-10)
            saved={'meta':meta,'v':v,'seconds':time.perf_counter()-t0,
                   'rows':[encode(k,a) for k,a in sorted(out.items(),key=lambda it:repr(it[0]))]}
            atomic_json(f,saved)
        nw=sum(abs(a)**2 for a in out.values());weight+=nw;nrows+=len(out);seen.update(out)
        gc=sum(abs(a)**2 for k,a in out.items() if 0 in k[0]);changed+=gc
        by_vertex.append({'v':v,'support':len(out),'norm2':nw,'graph_changed_norm2':gc})
        for k,a in sorted(out.items(),key=lambda it:-abs(it[1]))[:reverse_per_vertex]:
            if k==source:continue
            rp=folder/f'reverse_v{v}_target_{hashlib.sha256(repr(k).encode()).hexdigest()[:12]}.json'
            if rp.exists():
                obj=json.loads(rp.read_text())
                if obj['meta']!=meta or obj['target']!=encode(k,a):raise RuntimeError('stale reverse cache')
                b=complex(float.fromhex(obj['re']),float.fromhex(obj['im']))
            else:
                t0=time.perf_counter()
                rev=K5.prune_state(K5.apply_H_cached_state({k:1+0j},v,jmax2),1e-10)
                b=rev.get(source,0j)
                atomic_json(rp,{'meta':meta,'target':encode(k,a),'v':v,
                       're':float(b.real).hex(),'im':float(b.imag).hex(),
                       'support':len(rev),'seconds':time.perf_counter()-t0})
            reverse.append({'v':v,'conjugacy_error':abs(b-a.conjugate())})
        print(f'v={v} support={len(out)} norm2={nw:.12f}',flush=True)
    result={'provenance':meta,'scope':'one input, selected K5 H_v, finite wall; not physical projector',
       'completed_nodes':vertices,'num_direct_rows_with_duplicates':nrows,
       'unique_output_states':len(seen),'norm2':weight,'changed_graph_norm2':changed,
       'changed_graph_norm2_fraction':changed/weight if weight else 0.,
       'direct':by_vertex,'selected_reverse':reverse,
       'outside_input_shell_warning':'for source dimension 1, all amplitude may leave shell',
       'physical_projector_computed':False}
    atomic_json(folder/'results.json',result)
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--jmax2',type=int,default=5)
    p.add_argument('--source-index',type=int,default=0)
    p.add_argument('--vertices',type=int,nargs='+',default=list(range(5)))
    p.add_argument('--reverse-per-vertex',type=int,default=0)
    a=p.parse_args()
    r=run(a.output_dir,a.jmax2,a.source_index,a.vertices,a.reverse_per_vertex)
    print(json.dumps(r,indent=2))
if __name__=='__main__':main()
