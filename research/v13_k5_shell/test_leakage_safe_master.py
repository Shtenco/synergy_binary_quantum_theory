#!/usr/bin/env python3
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from bqg_k5_leakage_safe_master import assemble_columns

class LeakageSafeTests(unittest.TestCase):
    def test_false_null_from_compression(self):
        d={(0,"a"):{"outside":1.},(0,"b"):{}}
        r=assemble_columns(["a","b"],[0],lambda k,v:d[v,k],lambda k:k)
        self.assertEqual(r["spectrum"]["compressed_nullity_at_threshold"],2)
        self.assertEqual(r["spectrum"]["full_nullity_at_threshold"],1)
        self.assertEqual(r["spectrum"]["false_null_modes_eliminated"],1)
    def test_physical_cross_sector_as_finite_control_only(self):
        d={(0,"a"):{"outside":1.},(0,"b"):{"outside":-1.}}
        r=assemble_columns(["a","b"],[0],lambda k,v:d[v,k],lambda k:k)
        self.assertEqual(r["spectrum"]["full_nullity_at_threshold"],1)
        self.assertEqual(len(r["spectrum"]["finite_shell_cross_graph_blocks"]),2)
    def test_psd_leakage(self):
        d={(0,"a"):{"a":1.,"outside":1.},(0,"b"):{"b":2.}}
        r=assemble_columns(["a","b"],[0],lambda k,v:d[v,k],lambda k:k)
        self.assertAlmostEqual(r["total_column_norm2"],6.)
        self.assertAlmostEqual(r["discarded_norm2_in_old_compression"],1.)
        self.assertGreaterEqual(r["spectrum"]["omitted_psd_min_eigenvalue"],-1e-12)
    def test_direct_reverse_element(self):
        d={(0,"a"):{"outside":1j},(0,"outside"):{"a":-1j}}
        r=assemble_columns(["a"],[0],lambda k,v:d[v,k],lambda k:k,
            reverse_image=lambda k,v:d[v,k],reverse_limit=1)
        self.assertLess(r["reverse_pairs"][0]["direct_reversed_error"],1e-12)
    def test_large_does_not_claim_projector(self):
        r=assemble_columns(["a","b"],[0],lambda k,v:{},lambda k:k,dense_limit=1)
        self.assertFalse(r["spectrum"]["computed"])
    def test_vacuum_null(self):
        r=assemble_columns(["a"],[0],lambda k,v:{},lambda k:k)
        self.assertEqual(r["spectrum"]["full_nullity_at_threshold"],1)
        self.assertEqual(r["spectrum"]["projector_residual_norm"],0.)
if __name__=="__main__":unittest.main()
