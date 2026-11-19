import unittest
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import bqg_physical_history_decoherence_engine as e

class PhysicalOverlap(unittest.TestCase):
    def test_positive_physical_gram(self):
        eta=np.array([[2,.2j,0],[-.2j,1.5,.1],[0,.1,1]],complex)
        states=[np.array([1,0,1j]),np.array([.5j,1,1]),np.array([1,1j,.5])]
        r=e.physical_overlap_matrix(states,eta)
        self.assertTrue(r["positive_semidefinite"])
        self.assertGreaterEqual(r["minimum_eigenvalue"],-1e-10)
    def test_reject_indefinite_form(self):
        with self.assertRaises(ValueError):
            e.physical_overlap_matrix([np.array([1,0])],np.diag([1.,-1.]))

class InterferenceAndEnvironment(unittest.TestCase):
    def test_environment_suppresses_cross_term_both_weights_nonzero(self):
        k0=np.array([[1],[0]],complex); k1=np.array([[0],[1]],complex)
        F0=k0@k0.conj().T; F1=k0@k1.conj().T
        e0=np.array([[1,0],[0,0]],complex); e1=np.array([[0,0],[1,0]],complex)
        rho=np.full((2,2),.5,complex)
        _,D=e.decoherence_functional({"a":np.kron(F0,e0),"b":np.kron(F1,e1)},np.kron(rho,e0))
        self.assertAlmostEqual(D[0,0].real,.5)
        self.assertAlmostEqual(D[1,1].real,.5)
        self.assertAlmostEqual(abs(D[0,1]),0.)
    def test_without_environment_record_cross_term_survives(self):
        k0=np.array([[1],[0]],complex); k1=np.array([[0],[1]],complex)
        F0=k0@k0.conj().T; F1=k0@k1.conj().T
        _,D=e.decoherence_functional({"a":F0,"b":F1},np.full((2,2),.5,complex))
        self.assertAlmostEqual(D[0,1].real,.5)
    def test_support_and_coherent_depth_differ(self):
        r=e.history_sector_report([(1,"s",1+1j),(1,"s",-1-1j),(2,"s",.2j)])
        self.assertEqual(r["first_support_depth"],1)
        self.assertEqual(r["first_coherent_depth"],2)

class Refinement(unittest.TestCase):
    def test_three_channels_compared_per_refinement(self):
        levels={}
        for name,h in (("coarse",1.),("fine",.25)):
            v=2+.2*h
            levels[name]={c:[(.5,.5*v,1),(.5,.75*v,0),(1.,v,1),(1.,1.25*v,0),(1.5,1.5*v,1)]
                          for c in ("photon","matter","tensor")}
        r=e.refinement_channel_frontier_report(levels,1e-12)
        self.assertTrue(all(abs(x["cross_channel_spread"] or 0)<1e-10 for x in r["levels"]))
        self.assertEqual(r["status"],"FINITE_REFINEMENT_DIAGNOSTIC_ONLY")

if __name__ == "__main__":
    unittest.main(verbosity=2)
