"""V24 GitHub self-contained verifier tests. These are FINITE math checks."""
import math
import unittest
from fractions import Fraction as F
from verify_finite_gates import verify,source,inertia,complex_fraction,decimal_floor_sqrt,BRACKETS

class FiniteV24(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.r=verify()
  cls.s=source()
  cls.g=[[complex_fraction(z) for z in row] for row in cls.s['complex_gram_re_im']]
 def test_01_five_exact_intervals(self):
  self.assertEqual(len(self.r['W5_eigen_brackets']),5)
 def test_02_five_positive_pivots(self):
  self.assertTrue(all(z>0 for z in inertia(self.g,'0')[1]))
 def test_03_exact_inertia(self):
  for i,(lo,hi) in enumerate(BRACKETS):
   self.assertEqual(inertia(self.g,lo)[0],i)
   self.assertEqual(inertia(self.g,hi)[0],i+1)
 def test_04_strict_floor_sqrt(self):
  self.assertEqual(decimal_floor_sqrt(F(2)),'1.414213562373')
 def test_05_positive_k_rank(self):
  self.assertEqual([x['k'] for x in self.r['W5_kernel_nullity_Q_necessary_bounds']],[1,2,3,4,5])
 def test_06_w17_false_compressed_kernel(self):
  self.assertEqual(self.r['W17_compressed_zero_count'],15)
  self.assertEqual(self.r['W17_full_E0_kernel_dim'],0)
 def test_07_chiral_w43(self):
  self.assertEqual(self.r['W43_chiral_index_min_zero'],11)
  self.assertEqual(self.r['W43_chiral_zero_count'],11)
 def test_08_chiral_plusminus(self):
  self.assertLess(self.r['W43_nonzero_plus_minus_eigen_pairing_error'],1e-13)
 def test_09_hda_rank(self):
  self.assertEqual(self.r['HDA_ten_probe_rank'],10)
  self.assertEqual(self.r['HDA_six_distance_two_nonEE_required_rank'],6)
 def test_10_not_full_physics(self):
  self.assertFalse(self.r['physical_HDA_proved'])
  self.assertFalse(self.r['physical_continuum_proved'])
  self.assertFalse(self.r['interacting_TT_proved'])
if __name__=='__main__':unittest.main(verbosity=2)
