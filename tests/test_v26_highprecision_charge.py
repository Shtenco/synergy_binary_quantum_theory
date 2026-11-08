#!/usr/bin/env python3
"""Synthetic regression tests only: no physical SU(2) matrix elements here."""
from __future__ import annotations
import math,sys,tempfile,unittest
from pathlib import Path
from fractions import Fraction
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(ROOT/"research"/"v26"))
from v26_direct_branch_norm import squared_norm,exact_small_squared_norm,CertificationUnavailable
from v26_charge_scalar_diagnostics import load,diagnostics,minimum_consistent_uniform_eta,bound_absolute_charge_leak
from v26_hda_habitat import replay as habitat

def one(c,spin=(0,),v=(0,1,2)):
    return (spin,{i:[1+0j] for i in v},c)

class NumericRegression(unittest.TestCase):
    def test_01_exact_phase(self):
        self.assertEqual(mp.mpf(squared_norm([one(1+2j)],(0,1,2))["squared_norm"]),5)
    def test_02_orthogonal_sectors(self):
        a=[one(1,(0,)),one(3j,(1,))]
        self.assertEqual(mp.mpf(squared_norm(a,(0,1,2))["squared_norm"]),10)
    def test_03_cancellation(self):
        self.assertEqual(mp.mpf(squared_norm([one(1),one(-1)],(0,1,2))["squared_norm"]),0)
    def test_04_fraction_crosscheck(self):
        a=[one(1e8),one(1e-8),one(-1e8)]
        y=squared_norm(a,(0,1,2),digits=90)
        ex=exact_small_squared_norm(a,(0,1,2))
        with mp.workdps(95):
            self.assertTrue(mp.almosteq(mp.mpf(y["squared_norm"]),
                        mp.mpf(ex.numerator)/ex.denominator,rel_eps=mp.mpf("1e-70")))
    def test_05_no_fake_physical_claim(self):
        d=squared_norm([one(1)],(0,1,2))
        self.assertFalse(d["full_HDA_proved"])
        self.assertFalse(d["coherent_charge_projection_proved"])
    def test_06_tiny_amplitude_preserved(self):
        d=squared_norm([one(1e-200)],(0,1,2))
        self.assertGreater(mp.mpf(d["squared_norm"]),0)
    def test_07_nonfinite_rejected(self):
        with self.assertRaises(ValueError):squared_norm([one(float("nan"))],(0,1,2))
    def test_08_basis_limit_fail_closed(self):
        a=[((0,),{0:[1]*20,1:[1]*20,2:[1]*20},1)]
        with self.assertRaises(CertificationUnavailable):squared_norm(a,(0,1,2),max_group_dim=1000)
    def test_09_mixed_local_shape(self):
        a=[one(1),((0,),{0:[1,1],1:[1],2:[1]},1)]
        with self.assertRaises(ValueError):squared_norm(a,(0,1,2))
    def test_10_rational_work_limit(self):
        with self.assertRaises(CertificationUnavailable):
            exact_small_squared_norm([one(1)]*10,(0,1,2),max_terms=5)
    def test_11_actual_scalar_census(self):
        self.assertEqual(len(load()),12)
    def test_12_noncertified_three(self):
        d=diagnostics(load());self.assertEqual(d["noncertified"],[0,3,6])
    def test_13_minimum_error_consistency(self):
        self.assertAlmostEqual(minimum_consistent_uniform_eta(load()),3.0531133177191805e-16)
    def test_14_no_inconsistent_eta_zero(self):
        with self.assertRaisesRegex(ValueError,"NO_COMMON_ORTHOGONAL_PROJECTOR"):
            bound_absolute_charge_leak(load(),0)
    def test_15_only_conditional_charge_bound(self):
        d=load()
        self.assertGreater(bound_absolute_charge_leak(d,minimum_consistent_uniform_eta(d)),2.9e-7)
    def test_16_32_minimal_onehit_habitat(self):
        self.assertEqual(habitat()["minimal_combined_one_hit_habitat_dim_at_least"],32)
    def test_17_no_full_HDA(self):
        self.assertFalse(habitat()["full_HDA_proved"])

if __name__=="__main__":unittest.main()
