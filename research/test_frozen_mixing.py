"""Verify analytic predictions and unbiased sampling independently of plots."""
import unittest
import numpy as np
from frozen_mixing import estimates, moments


class FrozenMixingTests(unittest.TestCase):
    def test_exact_endpoint_moments(self):
        truth, _, vf, vd = moments(0, 2, 1)
        np.testing.assert_allclose(truth, [0, 1, np.exp(-.5)])
        np.testing.assert_allclose(vf[:2], [1, 2])
        np.testing.assert_allclose(vd, [0, 0, 0])
        truth, _, vf, vd = moments(1, 2, 1)
        np.testing.assert_allclose(truth[:2], [2, 5])
        np.testing.assert_allclose(vf[:2], [1, 18])
        np.testing.assert_allclose(vd[:2], [2, 20])

    def test_sampling_matches_exact_means_and_variances(self):
        truth, _, vf, vd = moments(.5, 2, 1)
        f, d, _, _ = estimates(12000, 64, 32, .5, 2, 1, [100, 200, 300], 256)
        for values, prediction in ((f, vf/64), (d, vd/32)):
            np.testing.assert_allclose(values.var(axis=0, ddof=1), prediction, rtol=.06)
            errors = np.abs(values.mean(axis=0)-truth)/np.sqrt(prediction/len(values))
            self.assertTrue(np.all(errors < 5))
        covariance = np.cov(f[:,0], d[:,0], ddof=1)[0,1]
        self.assertLess(abs(covariance), 5*np.sqrt(vf[0]/64*vd[0]/32/len(f)))

    def test_fixed_budget_cannot_beat_best_single_exact_variance(self):
        for epsilon in (.05, .25, .5, .9):
            _, _, vf, vd = moments(epsilon, 2, 1)
            a, b = vf/512, vd/256
            mixed = a*b/(a+b)
            best_single_budget = np.minimum(vf/1024, vd/512)
            self.assertTrue(np.all(mixed >= best_single_budget-1e-14))


if __name__ == "__main__":
    unittest.main()
