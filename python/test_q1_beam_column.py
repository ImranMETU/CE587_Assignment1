import unittest
import numpy as np
import q1_beam_column as q
class TestQ1(unittest.TestCase):
    def test_units(self):
        np.testing.assert_allclose(q.cm_to_m([2,.1]),[.02,.001])
        np.testing.assert_allclose(q.m_to_cm([.02,.001]),[2,.1])
        self.assertAlmostEqual(q.analytic()['mass_coefficient_kg_cm2'],14.702653618800232,11)
    def test_reference_formulas(self):
        x=q.cm_to_m([2,1/(2*np.pi)]);r=q.responses(x)
        self.assertAlmostEqual(r['area_m2'],.0002,14);self.assertAlmostEqual(r['inertia_m4'],4e-8,15)
        self.assertAlmostEqual(r['self_weight_N_m'],15.288,11)
        self.assertAlmostEqual(r['deflection_m'],.00191953125,14)
        self.assertAlmostEqual(r['stress_Pa']/1e6,250,10)
    def test_independent_and_cancellation(self):
        for x in [q.cm_to_m([2,.2]),q.cm_to_m([10,.5]),q.analytic()['representative_m']]:
            a=q.responses(x);b=q.independent_check(x)
            for k,v in a.items():self.assertAlmostEqual(v/b[k],1,12)
            self.assertAlmostEqual(b['deflection_m'],b['simplified_deflection_m'],14)
        self.assertAlmostEqual(q.deflection([.03,.002]),q.deflection([.03,.008]),14)
    def test_entire_optimal_segment(self):
        a=q.analytic()
        for R in np.linspace(*a['optimal_R_interval_cm'],1001):
            x=q.cm_to_m([R,a['product_cm2']/R]);self.assertTrue(q.feasible(x,tol=1e-12))
            self.assertAlmostEqual(q.mass(x),4.68,12)
        self.assertTrue(q.analytic()['deflection_radius_limit_cm']<2)
        self.assertLess(q.deflection([.02,.001]),q.DELTA_MAX)
    def test_activity_and_distinct_optima(self):
        right=q.analytic()['representative_m'];left=q.cm_to_m([2,1/(2*np.pi)])
        self.assertNotAlmostEqual(left[0],right[0],5)
        self.assertAlmostEqual(q.mass(left),q.mass(right),12)
        self.assertEqual({k for k,v in q.activity(right).items() if v=='active'},{'axial stress','R/t'})
        self.assertEqual({k for k,v in q.activity(left).items() if v=='active'},{'axial stress','minimum radius'})
    def test_invalid_and_violations(self):
        for x in [[0,.001],[.02,0],[np.nan,.002]]:self.assertFalse(q.feasible(x))
        self.assertFalse(q.feasible(q.cm_to_m([2,.1])))
        self.assertFalse(q.feasible(q.cm_to_m([3,.1])))
    def test_multistart_nonuniqueness(self):
        out=q.optimize();self.assertEqual(out['accepted_runs'],25)
        lo,hi=out['numerical_optimal_R_range_cm'];self.assertGreater(hi-lo,.1)
        self.assertTrue(any(r['initial_feasible'] for r in out['runs']))
        self.assertTrue(any(not r['initial_feasible'] for r in out['runs']))
if __name__=='__main__':unittest.main()
