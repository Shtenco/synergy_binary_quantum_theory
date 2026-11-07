#!/usr/bin/env python3
"""V21 pure-C1 structural regression tests; mock amplitudes are NOT physics."""
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from v21_stage_logic import all_32_specs,first_leg_census,operations,predecessor_ok
from v21_checkpoint_guard import validate_physical_stage
from pl_v21_pair_collector import compare_states,read_stage,read_v2
from pl_v20_stage_worker import save_charged,load_charged


class V21C1Tests(unittest.TestCase):
    def setUp(self):
        self.st={((1,)*32,(0,)*16,0,0,1,1):1+2j}
        self.one={'stage':1,'passed':True,'source_commit':'a'*40,'source_node':0,
                  'omit':0,'cycle':0,'mode':'forward','order':'p','aux':[0,0,1],
                  'Jmax2':7,'volume_backend':'prod','volume_backend_sha256':'b'*64,
                  'input_state_sha256':'c'*64,'operator_signature':['CV',4,1,0]}
        self.two=dict(self.one,stage=2,order='q',aux=[0,1,1],
                      operator_signature=['CK',2,1,1],post_second_scalar_filter=True)
        self.three=dict(self.two,stage=3,operator_signature=['CK',3,0,1])
    def test_01_census32(self):
        self.assertEqual(len(list(all_32_specs((2,3,4)))),32)
        self.assertEqual(first_leg_census((2,3,4)),(12,32))
    def test_02_original_forward(self):
        self.assertEqual(operations((2,3,4),'forward',(0,1,1)),
                         [('CV',4,1,0),('CK',3,1,1),('CK',2,0,1)])
    def test_03_original_adjoint(self):
        self.assertEqual(operations((2,3,4),'adjoint',(0,1,1)),
                         [('CK',2,1,0),('CK',3,1,1),('CV',4,0,1)])
    def test_04_safe_reuse_first(self):
        self.assertTrue(predecessor_ok(self.one,self.two,2,['CV',4,1,0]))
    def test_05_reject_wrong_first(self):
        with self.assertRaisesRegex(ValueError,'MISMATCH_FIRST_OPERATOR_SIGNATURE'):
            predecessor_ok(self.one,self.two,2,['CV',4,0,1])
    def test_06_strict_second_aux(self):
        with self.assertRaisesRegex(ValueError,'MISMATCH_aux'):
            predecessor_ok(self.two,dict(self.three,aux=[1,1,1]),3,None,['CK',2,1,1])
    def test_07_requires_scalar_filter(self):
        with self.assertRaisesRegex(ValueError,'STAGE2_FILTER_NOT_APPLIED'):
            predecessor_ok(dict(self.two,post_second_scalar_filter=False),
                           self.three,3,None,['CK',2,1,1])
    def test_08_reject_backend_mismatch(self):
        with self.assertRaisesRegex(ValueError,'MISMATCH_volume_backend_sha256'):
            predecessor_ok(self.one,dict(self.two,volume_backend_sha256='f'*64),2,['CV',4,1,0])
    def test_09_low_physical_leak(self):
        self.assertTrue(validate_physical_stage({'CV_complete_basis_leakage':1e-14},self.st,1)['stage_physical_basis_pass'])
    def test_10_high_physical_leak(self):
        with self.assertRaisesRegex(ValueError,'PHYSICAL_BASIS_LEAKAGE'):
            validate_physical_stage({'internal_volume_sector_leakage':1.},self.st,2)
    def test_11_charge_leak_separate(self):
        d={'outer_complete_basis_leakage':0.,'complete_charge_basis_leakage':1.}
        x=validate_physical_stage(d,self.st,2)
        self.assertEqual(x['charge_leakage_max_unprojected'],1.)
        self.assertFalse(x['charge_covariance_verified'])
    def test_12_missing_physical_diagnostic(self):
        with self.assertRaisesRegex(ValueError,'MISSING_PHYSICAL_LEAKAGE'):
            validate_physical_stage({'complete_charge_basis_leakage':0.},self.st,3)
    def test_13_empty_not_physical_proof(self):
        x=validate_physical_stage({'CV_complete_basis_leakage':0.},{},1)
        self.assertTrue(x['stage_zero_support'])
        self.assertFalse(x['charge_covariance_verified'])
    def test_14_complex_phase_matters(self):
        k=next(iter(self.st))
        self.assertGreater(compare_states({k:1+0j},{k:0+1j})['relative_delta'],1.4)
    def test_15_checkpoint_integrity(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.npz';h=save_charged(p,self.st)
            self.assertEqual(load_charged(p,h),self.st)
            with self.assertRaises(ValueError):load_charged(p,'0'*64)
    def test_16_collector_rejects_unguarded_stage(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'stage';h=save_charged(p.with_suffix('.npz'),self.st)
            p.with_suffix('.json').write_text(json.dumps(
                dict(self.one,passed=True,status='REAL_ORIGINAL_STAGE_ONLY',
                     npz_sha256=h,state_support=1)))
            with self.assertRaisesRegex(ValueError,'MISSING_V21_ACCEPTANCE_GUARD'):
                read_stage(p)
    def test_17_read_original_v2(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'v2.npz'
            np.savez_compressed(p,spins=np.array([[1]*32]),Ks=np.array([[0]*16]),
                                amp=np.array([1+2j]))
            sha=hashlib.sha256(p.read_bytes()).hexdigest()
            self.assertEqual(list(read_v2({'npz':str(p),'sha256':sha}).values()),[1+2j])
    def test_18_reject_v2_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'v2.npz';np.savez_compressed(p,spins=np.empty((0,32)),
                               Ks=np.empty((0,16)),amp=np.empty((0,),dtype=complex))
            with self.assertRaisesRegex(ValueError,'INDEPENDENT_V2_SHA256_MISMATCH'):
                read_v2({'npz':str(p),'sha256':'0'*64})


if __name__=='__main__':
    unittest.main()
