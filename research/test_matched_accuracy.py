import unittest
import numpy as np
from matched_accuracy import choose,initial,metrics,reference_estimates,scales


class MatchedAccuracyTests(unittest.TestCase):
    def test_reference_control_variate_removes_constant_initial_noise(self):
        truth=initial(.05)
        noise=np.array([[.1,.2,.3],[-.1,-.2,-.3]])[:,None,:]
        samples=np.broadcast_to(truth+noise,(2,3,3)).copy()
        reference=reference_estimates({'estimates':{'full':samples}},.05)
        np.testing.assert_allclose(reference,np.broadcast_to(truth,(2,3,3)),atol=1e-14)

    def test_selection_uses_error_upper_bound_and_measured_runtime(self):
        rows=[{'label':'uncertain','metrics':{'full':{'joint_relative_rmse_ci95':[.1,.5],'mean_compute_seconds':.01}}},
              {'label':'slow','metrics':{'full':{'joint_relative_rmse_ci95':[.1,.2],'mean_compute_seconds':.3}}},
              {'label':'fast','metrics':{'full':{'joint_relative_rmse_ci95':[.1,.2],'mean_compute_seconds':.1}}}]
        self.assertEqual(choose(rows,'full',.25)['label'],'fast')
        self.assertIsNone(choose(rows,'full',.05))

    def test_normalized_trajectory_error_matches_known_offset(self):
        truth=np.broadcast_to(initial(.05),(8,3,3)).copy()
        estimates=np.broadcast_to(initial(.05)+scales(.05),(8,3,3)).copy()
        data={'settings':{'mode':'pic'},'estimates':{'full':estimates,'pic_cv':truth},'costs':np.ones(8)*.1}
        result=metrics(data,truth,.05,100,bootstraps=20)
        self.assertAlmostEqual(result['full']['joint_relative_rmse'],1)
        self.assertAlmostEqual(result['pic_cv']['joint_relative_rmse'],0)


if __name__=='__main__':unittest.main()
