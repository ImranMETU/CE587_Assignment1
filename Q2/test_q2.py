import unittest
import numpy as np
import q2_cantilever as q
class TestQ2(unittest.TestCase):
    def test_solid(self):
        r=.1;self.assertAlmostEqual(q.geometry([r,0])[0],np.pi*r*r)
        self.assertAlmostEqual(q.geometry([r,0])[1],np.pi*r**4/4)
    def test_invalid(self):
        for x in [[0,0],[.2,.2],[.1,.2],[.2,-.1],[np.nan,0]]:
            self.assertFalse(q.feasible(x))
            with self.assertRaises(ValueError):q.geometry(x)
    def test_independent(self):
        x=q.analytical_optimum();r=q.independent_check(x)
        self.assertAlmostEqual(q.mass(x),r['mass'],8)
        self.assertAlmostEqual(q.bending_stress(x)/r['bending'],1,12)
        self.assertAlmostEqual(q.shear_stress(x)/r['shear'],1,12)
        self.assertTrue(q.feasible(x))
        self.assertAlmostEqual(r['mass'],334.447684923,6)
    def test_governing_and_global(self):
        C=4*q.P*q.L/(np.pi*q.SIGMA)
        ro=np.linspace(C**(1/3),.4,1000)
        ri=np.maximum(ro**4-C*ro,0)**.25
        masses=np.array([q.mass(x) for x in zip(ro,ri)])
        self.assertTrue(np.all(np.diff(masses)<0))
        self.assertTrue(all(q.shear_stress(x)<=6.6e6*(1+1e-12) for x in zip(ro,ri)))
    def test_multistart_and_guard(self):
        out=q.optimize();self.assertEqual(out['accepted_runs'],20)
        self.assertEqual(set(out['active']),{'bending','outer maximum'})
        for guard in [1e-10,1e-6]:
            r=q.optimize(guard);self.assertAlmostEqual(r['mass_kg'],out['mass_kg'],4)
    def test_violations(self):
        self.assertFalse(q.feasible([.4,.3999]))
        self.assertFalse(q.feasible([.41,.3]))
if __name__=='__main__':unittest.main()
