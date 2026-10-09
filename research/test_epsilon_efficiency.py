import unittest
import numpy as np
from epsilon_efficiency import combine, ratio_interval


class EfficiencyTests(unittest.TestCase):
    def test_ratio_orientation_and_batch_bootstrap(self):
        pic=[np.full(8,2.) for _ in range(3)]
        hdp=[np.full(8,.5) for _ in range(3)]
        np.testing.assert_allclose(ratio_interval(pic,hdp,17,100),[.25,.25])

    def test_batches_preserve_whole_trajectories(self):
        batches=[dict(settings={'mode':'pic'},costs=np.full(2,i),
                      estimates={'full':np.full((2,4,3),i)}) for i in (1,2,3)]
        result=combine(batches)
        self.assertEqual(result['estimates']['full'].shape,(6,4,3))
        np.testing.assert_array_equal(result['estimates']['full'][:,0,0],[1,1,2,2,3,3])
        np.testing.assert_array_equal(result['costs'],[1,1,2,2,3,3])


if __name__=='__main__':
    unittest.main()
