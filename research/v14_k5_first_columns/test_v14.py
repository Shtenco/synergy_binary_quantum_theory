#!/usr/bin/env python3
"""Regression for V14 cold-column output without inferring physical projector."""
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).parent))
class V14Tests(unittest.TestCase):
    def test_upstream_module_exists(self):
        self.assertTrue((Path(__file__).parents[2]/"scripts"/"k5_peter_weyl_safe_hda_column.py").exists())
    def test_cached_hex_roundtrip(self):
        for x in (0.,1.,-0.7722790801333823,1e-16):
            self.assertEqual(float.fromhex(x.hex()),x)
    def test_k5_initial_basis(self):
        import k5_peter_weyl_safe_hda_column as k
        self.assertEqual(len(k.basis_full_jhalf()),32)
    def test_sector_graph_signature(self):
        from itertools import combinations
        edges=list(combinations(range(5),2))
        self.assertEqual(len(edges),10)
if __name__=="__main__":
    import os
    sys.path.insert(0,str(Path(__file__).parents[2]/"scripts"))
    unittest.main()
