#!/usr/bin/env python3
"""V20 regression tests: printed-data math + mock-only staged bookkeeping."""
from __future__ import annotations
import copy
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"research"/"v20"))
sys.path.insert(0,str(ROOT/"scripts"))
from v20_rank_thresholds import source_read, exact_printed_gershgorin, replay
from pl_v20_stage_worker import (ordered_operations,scalar_filter,do_leg,save_charged,
                                 load_charged,validate_predecessor)

SOURCE=ROOT/"research"/"v20"/"SOURCE_SNAPSHOT.json"


class V20Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.src=source_read(SOURCE)
        cls.out=replay(SOURCE)

    def test_01_provenance(self):
        self.assertEqual(self.src["source_git_blob_sha"],
                         "0231c54d2e9a7fc9401995c999d3989f3f1f7648")
    def test_02_exact_printed_positive(self):
        self.assertGreater(exact_printed_gershgorin(self.src),30)
    def test_03_matrix_hermitian(self):
        changed=copy.deepcopy(self.src)
        changed["complex_gram_re_im"][0][1][0]+=1
        with self.assertRaises(ValueError):
            exact_printed_gershgorin(changed)
    def test_04_reference_spectrum(self):
        np.testing.assert_allclose(
            self.out["K5_W5"]["eigenvalues"],
            [float(z) for z in self.src["published_eigvals"]],atol=1e-11)
    def test_05_rank_k_monotonic(self):
        x=self.out["K5_W5"]["kernel_nullity_compensator_thresholds"]
        self.assertEqual([r["k"] for r in x],list(range(1,6)))
        self.assertTrue(all(x[i]["necessary_Q_Frobenius_norm"]<
                            x[i+1]["necessary_Q_Frobenius_norm"]
                            for i in range(4)))
    def test_06_reject_corrupt_source_sha(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"source.json"
            f.write_bytes(SOURCE.read_bytes()+b" ")
            with self.assertRaisesRegex(ValueError,"SOURCE_SHA256_MISMATCH"):
                source_read(f)
    def test_07_no_physics_claim(self):
        self.assertFalse(self.out["new_physical_Lorentzian_amplitudes"])
        self.assertFalse(self.out["full_HDA_proved"])
    def test_08_rank_obstruction(self):
        x=self.out["PL16_ten_EE_probes"]["rank_limited_compensator_residual_floors"]
        self.assertEqual(len(x),10)
        self.assertAlmostEqual(
            x[-1]["necessary_HDA_residual_operator_lower"],4/9*1.324189319007516)
    def test_09_forward_order(self):
        self.assertEqual(
            ordered_operations((10,11,12),"forward",(1,0,1)),
            [("CV",12,1,1),("CK",11,0,1),("CK",10,1,0)])
    def test_10_adjoint_order(self):
        self.assertEqual(
            ordered_operations((10,11,12),"adjoint",(1,0,1)),
            [("CK",10,0,1),("CK",11,1,0),("CV",12,1,1)])
    def test_11_exact_filter(self):
        src={((),(),x,0,0,0):1j for x in range(4)}
        self.assertEqual({q[2] for q in scalar_filter(src)},{0,2})
    def test_12_mock_operations_not_science(self):
        seen=[]
        def cv(state,src,t,l,r,J):
            seen.append(("CV",t,l,r,J))
            return {key:z*(2+1j) for key,z in state.items()},0.0
        def ck(state,src,t,l,r,J):
            seen.append(("CK",t,l,r,J))
            return {key:z*(1-3j) for key,z in state.items()},{"leakage":0.0}
        st={((1,)*32,(0,)*16,0,0,0,0):1+0j}
        for i,leg in enumerate(ordered_operations((2,3,4),"forward",(0,1,0)),1):
            st,_=do_leg(st,leg,0,cv,ck)
            if i==2:st=scalar_filter(st)
        self.assertEqual(seen,[("CV",4,0,0,7),("CK",3,1,0,7),("CK",2,0,1,7)])
    def test_13_complex_checkpoint(self):
        state={((1,)*32,(0,)*16,0,0,1,1):3+5j,
               ((2,)*32,(2,)*16,2,-1,0,0):1e-20+2j}
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"s.npz"
            digest=save_charged(f,state)
            self.assertEqual(load_charged(f,digest),state)
            with self.assertRaises(ValueError):
                load_charged(f,"f"*64)
    def test_14_mixed_backend_refused(self):
        prev={"source_commit":"a"*40,"source_node":0,"omit":0,"cycle":0,
              "mode":"forward","order":"p","aux":[0,1,0],"Jmax2":7,
              "volume_backend":"original","npz_sha256":"b"*64,
              "stage":2,"passed":True}
        validate_predecessor(prev,dict(prev,stage=3),3)
        with self.assertRaisesRegex(ValueError,"PREDECESSOR_MISMATCH_volume_backend"):
            validate_predecessor(prev,dict(prev,stage=3,volume_backend="other"),3)
    def test_15_nonfinite_rejected(self):
        key=((1,)*32,(0,)*16,0,0,0,0)
        def bad(*args):return {key:complex(float("nan"),0)},0.0
        with self.assertRaisesRegex(ValueError,"NONFINITE_AMPLITUDE"):
            do_leg({key:1+0j},("CV",0,0,0),0,bad,None)


if __name__=="__main__":
    unittest.main(verbosity=2)
