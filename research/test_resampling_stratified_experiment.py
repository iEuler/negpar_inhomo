import unittest
import numpy as np
from test_resampling_domain_experiment import domain_fixture
from resampling_stratified_experiment import analyze_stratified

class StratifiedTests(unittest.TestCase):
    def fixture(self):
        old=domain_fixture()
        data=np.zeros(len(old),dtype=old.dtype.descr+[('stratified','f8')])
        for n in old.dtype.names: data[n]=old[n]
        data['stratified']=data['aligned']; data['aligned']=data['cutoff']>0
        return data
    def test_paired_rmse_and_shared_full_control(self):
        s=analyze_stratified(self.fixture()[::-1],8,2)
        for method in ('independent','stratified'): self.assertEqual(len(s[method]['rows']),48)
        for r in s['allocation_pairs']:
            for v in r['rmse_ratios'].values():
                self.assertAlmostEqual(v['ratio'],.5)
                self.assertAlmostEqual(v['ci95'][0],.5)
                self.assertAlmostEqual(v['ci95'][1],.5)
    def test_invalid_domain_and_source_rejected(self):
        d=self.fixture(); d['aligned'][-1]=0
        with self.assertRaisesRegex(RuntimeError,'aligned core'): analyze_stratified(d,8,2)
        d=self.fixture(); d['source_0'][-1]=1
        with self.assertRaisesRegex(RuntimeError,'Source mismatch'): analyze_stratified(d,8,2)

if __name__=='__main__': unittest.main()
