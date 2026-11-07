#!/usr/bin/env python3
"""Synthetic ONLY: V22 stage-link, coherent-collection and provenance safeguards."""
from __future__ import annotations
import hashlib,json,sys,tempfile,types,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from pl_v22_stage_worker import save_charged
import pl_v22_pair_collector as C
from v21_stage_logic import operations,all_32_specs

SP=(1,)*32;CH=(SP,(-1,)+(0,)*15,0,0,0,0);G=(SP,(0,)*16)
class D:
    def __init__(self,*args):self.neighbor={(0,i):i+1 for i in range(4)}
def pair(*args):return (0,1,2),(1,0,2),0,1
def spec(D,source,index):return (0,(),(),(),1 if index==0 else -1)
def projector(state,source,tol):
    return {G:sum(state.values())},float(sum(abs(z)**2 for z in state.values())),0.
def add(out,inp,scale=1,tol=1e-10):
    for k,z in inp.items():out[k]=out.get(k,0j)+scale*z

class V22C1Integration(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.base=Path(self.tmp.name)
        modules={
            'pl_dual_complex':types.SimpleNamespace(DualComplex=D,seed_16cell_boundary=lambda:None),
            'pl_lorentzian_physical_pair_worker_v4':types.SimpleNamespace(pair_indices=pair,add=add),
            'pl_lorentzian_triple_worker':types.SimpleNamespace(ordered_spec=spec,add=add,TOL=1e-10,max_spin=lambda s:0),
            'peter_weyl_lorentzian_gauss_action_gate':types.SimpleNamespace(project_scalar_gauss=projector)
        }
        self.patch=patch.dict(sys.modules,modules);self.patch.start()
        self.manifest={'source_commit':'a'*40,'source_node':0,'omit':0,'cycle':0,
                       'Jmax2':7,'volume_backend_sha256':'b'*64,
                       'input_state_sha256':'c'*64,'paths':[],'independent_v2':{}}
        for mode in ('forward','adjoint'):
            for order in ('p','q'):
                value_sum=0j
                for i in range(2):
                    for j in range(2):
                        for k in range(2):
                            aux=[i,j,k]
                            op=operations((1,2,3) if order=='p' else (2,1,3),mode,aux)
                            val=complex(1+i+j*.4+k*.2,.1+i*.2)
                            if mode=='adjoint':val*=1j
                            if order=='q':val*=.5
                            value_sum+=val
                            name=f'{mode}_{order}_{i}{j}{k}';prev=None;entry={}
                            for stage in (1,2,3):
                                prefix=self.base/(name+str(stage))
                                npz_sha=save_charged(prefix.with_suffix('.npz'),{CH:val if stage==3 else 1+0j})
                                m=dict(self.manifest,stage=stage,mode=mode,order=order,aux=aux,
                                       volume_backend='prod',passed=True,
                                       status='REAL_ORIGINAL_STAGE_ONLY',npz_sha256=npz_sha,
                                       state_support=1,operator_signature=op[stage-1],
                                       post_second_scalar_filter=stage>=2,
                                       guard={'stage_physical_basis_pass':True})
                                if prev:
                                    m['predecessor_npz_sha256']=hashlib.sha256(prev.with_suffix('.npz').read_bytes()).hexdigest()
                                    m['predecessor_metadata_sha256']=hashlib.sha256(prev.with_suffix('.json').read_bytes()).hexdigest()
                                prefix.with_suffix('.json').write_text(json.dumps(m))
                                entry['stage'+str(stage)]=str(prefix);prev=prefix
                            self.manifest['paths'].append(entry)
                label=('f' if mode=='forward' else 'a')+('p' if order=='p' else 'q')
                p=self.base/(label+'_v2.npz')
                np.savez_compressed(p,spins=np.array([SP]),Ks=np.array([(0,)*16]),amp=np.array([value_sum]))
                m={'passed':True,'mode':mode,'index':0 if order=='p' else 1,'source_node':0,
                   'Jmax':3.5,'checks':{'ok':True},
                   'producer_kind':'v22_independent_original_v2_worker','V4_staged_code_invoked':False,
                   'source_commit':'a'*40,'volume_backend_sha256':'b'*64,
                   'npz_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
                q=self.base/(label+'_meta.json');q.write_text(json.dumps(m))
                self.manifest['independent_v2'][label]={
                    'npz':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
                    'metadata_json':str(q),'metadata_sha256':hashlib.sha256(q.read_bytes()).hexdigest()}
    def tearDown(self):self.patch.stop();self.tmp.cleanup()
    def test_complete_32_path_mock(self):
        r,s=C.collect(self.manifest)
        self.assertEqual(r['paths_completed'],32)
        self.assertFalse(r['HDA_proved'])
        self.assertLess(max(v['relative_delta'] for v in r['comparisons'].values()),1e-14)
    def test_missing_path_fails(self):
        self.manifest['paths'].pop()
        with self.assertRaisesRegex(ValueError,'V4_32_PATHS_NOT_COMPLETE'):C.collect(self.manifest)
    def test_wrong_v2_commit_fails(self):
        e=self.manifest['independent_v2']['fp'];p=Path(e['metadata_json'])
        m=json.loads(p.read_text());m['source_commit']='f'*40;p.write_text(json.dumps(m))
        e['metadata_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
        with self.assertRaisesRegex(ValueError,'V2_COMMIT_DOES_NOT_MATCH_V4'):C.collect(self.manifest)
    def test_broken_stage_link_fails(self):
        prefix=Path(self.manifest['paths'][0]['stage3']);p=prefix.with_suffix('.json')
        m=json.loads(p.read_text());m['predecessor_npz_sha256']='f'*64;p.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError,'BROKEN_PREDECESSOR_NPZ_CHAIN'):C.collect(self.manifest)
    def test_subthreshold_auxiliary_amplitudes_must_cohere(self):
        # Eight terms individually below 1e-10 can add to a nonzero amplitude.
        # This fixture tests algebra, NOT original SU2 physics.
        parts=[6e-11+0j]*8
        sequential={}
        for z in parts:
            a=sequential.get('key',0j)+z
            if abs(a)>1e-10:sequential['key']=a
            elif 'key' in sequential:del sequential['key']
        self.assertEqual(sequential,{})
        self.assertAlmostEqual(sum(parts).real,4.8e-10,places=20)
    def test_12_first_signatures(self):
        specs=list(all_32_specs((1,2,3)))
        self.assertEqual(len(specs),32)
        self.assertEqual(len({tuple(x['operations'][0]) for x in specs}),12)

if __name__=='__main__':unittest.main()
