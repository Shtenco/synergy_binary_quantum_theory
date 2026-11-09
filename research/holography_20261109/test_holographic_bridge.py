#!/usr/bin/env python3
"""Regression tests for finite holography-inspired bridges (not BQG physics)."""
import sys
import unittest
from pathlib import Path
from math import comb,log
from itertools import product,combinations
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_holographic_bridge import run

class FiniteHolographyBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out=run()
    def test_wall_microstates(self):
        self.assertEqual(comb(8,4),70)
        self.assertAlmostEqual(self.out['wall_missing_information_bits'],log(70,2))
    def test_codewords_are_unique(self):
        self.assertEqual(len(set(self.out['boundary_parity_codewords'])),4)
    def test_two_classical_flags(self):
        self.assertTrue(all(self.out['any_two_classical_flags_recover_label'].values()))
    def test_isometry(self):
        self.assertEqual(self.out['boundary_encoding_isometry_rank'],4)
    def test_operator_reconstruction_global(self):
        self.assertTrue(self.out['operator_reconstruction_on_global_code_subspace'])
    def test_not_erasure_correcting_quantum_code(self):
        self.assertFalse(self.out['quantum_erasure_correction_for_any_two_qubit_logical_state'])
    def test_walsh_dot(self):
        self.assertEqual(self.out['walsh_pairwise_dot'],'-1/3')
    def test_q4_graph(self):
        self.assertEqual((self.out['Q4_vertices'],self.out['Q4_edges']),(16,32))
    def test_q4_cut(self):
        self.assertEqual(self.out['Q4_coordinate_cut_edges'],8)
    def test_q4_rank_and_entropy(self):
        self.assertEqual(self.out['Q4_schmidt_rank_bound'],256)
        self.assertEqual(self.out['Q4_von_neumann_entropy_bound_bits'],8)
    def test_witness_for_saturation(self):
        self.assertTrue(self.out['Q4_bell_pair_saturation_pure_state'])
    def test_no_false_physical_claims(self):
        for key in ('original_BQG_S_U_2_amplitudes_computed',
                    'physical_Bekenstein_Hawking_law_derived_from_BQG',
                    'AdS_CFT_duality_established_for_BQG',
                    'physical_quantum_gravity_proved'):
            self.assertFalse(self.out[key],key)

if __name__=='__main__':unittest.main(verbosity=2)
