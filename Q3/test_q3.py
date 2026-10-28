import unittest
import numpy as np
import q3_tripod as q
class TestQ3(unittest.TestCase):
    def test_conversions(self):
        self.assertEqual(q.E,75e9/1e4);self.assertEqual(q.SIGMA,150e6/1e4)
        self.assertEqual(q.RHO,2800/1e6)
    def test_independent_geometry_equilibrium_and_units(self):
        for x in [[50,4],[125,6],[200,10],q.analytical_optimum()]:
            r=q.responses(x);s=q.independent_check(x)
            for k,v in r.items():self.assertAlmostEqual(v/s[k],1,12)
            self.assertTrue(np.allclose(s['base_sides_m'],1.2))
            self.assertLess(np.linalg.norm(s['equilibrium_residual_N']),1e-9)
    def test_safety_factor_once(self):
        x=q.analytical_optimum()
        self.assertAlmostEqual(q.critical_load(x)/q.axial_force(x[0]),2,11)
    def test_boundary_and_monotonicity(self):
        H=np.linspace(50,200,2000);db=q.diameter_buckling(H);ds=q.diameter_compression(H)
        self.assertTrue(np.all(db>ds));self.assertTrue(np.all((db>.5)&(db<10)))
        self.assertTrue(np.all(np.diff(q.mass([H,db]))>0))
        self.assertAlmostEqual(q.compressive_stress(q.analytical_optimum())/q.SIGMA,.24755322936,10)
    def test_invalid_and_infeasible(self):
        for x in [[0,2],[50,0],[np.nan,3]]:self.assertFalse(q.feasible(x))
        self.assertFalse(q.feasible([50,2]));self.assertFalse(q.feasible([49,5]))
    def test_multistart(self):
        r=q.optimize();self.assertEqual(r['accepted_runs'],25)
        self.assertEqual(set(r['active']),{'Euler buckling','minimum height'})
        self.assertAlmostEqual(r['responses']['mass_kg'],6.6054480656,8)
        self.assertLess(abs(r['reduced_mass_kg']-r['responses']['mass_kg']),1e-6)
if __name__=='__main__':unittest.main()
