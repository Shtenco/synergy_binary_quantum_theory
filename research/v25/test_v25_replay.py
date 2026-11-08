import unittest
from check_hda_geometry import run as hda
from check_charge_summary import run as charge
class Tests(unittest.TestCase):
    def test_one(self):self.assertEqual(hda()["PL16_W17"]["offhabitat_rank_at_least"],15)
    def test_two(self):self.assertAlmostEqual(hda()["PL16_W17"]["offhabitat_scaled_Frobenius_lower"],5.184991899890351)
    def test_three(self):self.assertAlmostEqual(hda()["PL16_HDA_ten_probes"]["if_all_nonEE_and_independent_D_compensation_inside_W17_residual_F_lower"],3.3815397373424494)
    def test_four(self):self.assertEqual(charge()["unstable_cancellation_orientations"],[0,3,6])
    def test_five(self):self.assertEqual(charge()["inexact_zero_projection_violations"],[0,3])
    def test_six(self):self.assertFalse(hda()["full_Lorentzian_HDA_proved"])
if __name__=="__main__":unittest.main()
