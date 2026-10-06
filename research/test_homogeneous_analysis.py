import unittest
import numpy as np
from homogeneous_experiment import fit


class HomogeneousAnalysisTests(unittest.TestCase):
    def test_covariance_aware_weight_accounts_for_shared_background(self):
        rng=np.random.default_rng(123)
        common=rng.normal(size=(50000,1,1))
        full=common+rng.normal(size=common.shape)
        signed=common+2*rng.normal(size=common.shape)
        weight,independent,vf,vd,cov=fit(dict(full=full,signed=signed))
        # V_full=2, V_signed=5, C=1: full weight=4/5.
        self.assertAlmostEqual(float(weight[0,0]),.8,delta=.015)
        self.assertAlmostEqual(float(independent[0,0]),5/7,delta=.015)
        mixed=weight*full+(1-weight)*signed
        self.assertLess(float(mixed.var()),float(full.var()))

    def test_identical_constant_estimators_have_safe_weight(self):
        values=np.ones((100,2,3))
        weight,independent,*_=fit(dict(full=values,signed=values))
        np.testing.assert_array_equal(weight,0)
        np.testing.assert_array_equal(independent,0)


if __name__=='__main__':unittest.main()
