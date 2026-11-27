#!/usr/bin/env python3
"""Independent synthetic regression for the V18 Schur algebra.
This does NOT reproduce physical K5 numerics; use published NPZ inputs for those.
"""
import importlib.util, json, tempfile, unittest
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("v18_gate", HERE/"v18_schur_audit.py")
v18=importlib.util.module_from_spec(spec);spec.loader.exec_module(v18)

class V18SchurRegression(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        init=[((1,)*10,tuple((2 if i & (1<<j) else 0) for j in range(5))) for i in range(32)]
        extra=[
            ((0,1,1,1,1,1,1,1,1,1),(0,0,0,0,0)),
            ((0,0,1,1,1,1,1,1,1,1),(0,0,0,0,0)),
            ((1,0,1,1,1,1,1,1,1,1),(0,0,0,0,0)),
            ((1,1,1,1,1,1,1,1,1,2),(0,0,0,0,0)),
        ]
        # extras positions 32:odd, 33:even, 34:odd, 35:odd
        self.keys=np.asarray([repr(x) for x in init+extra])
        M=np.zeros((36,36),complex)
        M[:32,:32]=10*np.eye(32)
        M[33,33]=9.5
        M[0,33]=0.4
        M[33,0]=0.4
        for i,z in zip((32,34,35),(7.,8.,12.)):M[i,i]=z
        self.M=M
        self.v16=self.root/"v16.npz";self.v17=self.root/"v17.npz";self.status=self.root/"status.json"
        np.savez_compressed(self.v16,M=M[:32,:32])
        np.savez_compressed(self.v17,M=M,states=self.keys)
        self.status.write_text(json.dumps({"backend_git_blob_sha":v18.MIRROR_SHA,"canonical_sha_confirmed":False}))
    def tearDown(self):self.temp.cleanup()
    def test_exact_schur(self):
        r=v18.audit(self.v16,self.v17,self.status)
        self.assertEqual(r["input_dimensions"]["even"],33)
        self.assertAlmostEqual(r["coupling_norm"],0.4)
        self.assertLess(r["feshbach_root_error"],1e-12)
        self.assertEqual(r["numeric_nullity"],0)
    def test_wrong_backend_refused(self):
        self.status.write_text(json.dumps({"backend_git_blob_sha":"wrong"}))
        with self.assertRaises(ValueError):v18.audit(self.v16,self.v17,self.status)
    def test_cross_parity_refused(self):
        M=self.M.copy();M[0,32]=M[32,0]=0.1
        np.savez_compressed(self.v17,M=M,states=self.keys)
        with self.assertRaises(ValueError):v18.audit(self.v16,self.v17,self.status)
    def test_inconsistent_base_refused(self):
        A=self.M[:32,:32].copy();A[1,1]+=0.1
        np.savez_compressed(self.v16,M=A)
        with self.assertRaises(ValueError):v18.audit(self.v16,self.v17,self.status)
    def test_nonhermitian_refused(self):
        M=self.M.copy();M[1,2]=1j
        np.savez_compressed(self.v17,M=M,states=self.keys)
        with self.assertRaises(ValueError):v18.audit(self.v16,self.v17,self.status)
    def test_no_physical_claim(self):
        r=v18.audit(self.v16,self.v17,self.status)
        self.assertFalse(r["canonical_backend_verified"])
        self.assertFalse(r["full_physical_projector"])
if __name__=="__main__":unittest.main(verbosity=2)
