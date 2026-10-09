import unittest
import numpy as np
from test_resampling_domain_experiment import domain_fixture
from resampling_correction_experiment import analyze_corrected

class CorrectionTests(unittest.TestCase):
    def fixture(self):
        old=domain_fixture()
        extra=[('corrected','f8'),('stratified','f8')]+[(n,'f8') for n in ('correction_status','correction_iterations','correction_removed','correction_rms','correction_max','correction_residual')]
        extra += [(f'{prefix}_{j}','f8') for j in range(7) for prefix in ('call_target','correction_moment')]
        data=np.zeros(len(old),dtype=old.dtype.descr+extra)
        for n in old.dtype.names: data[n]=old[n]
        data['corrected']=data['aligned']; data['aligned']=data['cutoff']>0
        data['stratified']=data['cutoff']>0; data['correction_status'][data['corrected']==0]=-1
        return data
    def test_pairs_and_feasibility(self):
        s=analyze_corrected(self.fixture()[::-1],8,2)
        for method in ('uncorrected','corrected'): self.assertEqual(len(s[method]['rows']),48)
        for r in s['correction_pairs']:
            for v in r['rmse_ratios'].values():
                self.assertAlmostEqual(v['ratio'],.5)
                self.assertAlmostEqual(v['ci95'][0],.5)
                self.assertAlmostEqual(v['ci95'][1],.5)
        row=next(r for r in s['corrected']['rows'] if r['cutoff']==3)
        self.assertEqual(row['correction_success_fraction'],1)
    def test_invalid_success_constraint_rejected(self):
        d=self.fixture(); d['correction_moment_6'][-1]=1
        with self.assertRaisesRegex(RuntimeError,'violates measured moment'): analyze_corrected(d,8,2)
    def test_failed_correction_is_retained_in_errors(self):
        d=self.fixture(); d['correction_status'][d['corrected']==1]=6
        d['correction_moment_6'][d['corrected']==1]=1
        s=analyze_corrected(d,8,2)
        self.assertTrue(all(r['correction_success_fraction']==0 for r in s['corrected']['rows']))
        self.assertTrue(all(abs(r['rmse_ratios']['fourier']['ratio']-.5)<1e-12 for r in s['correction_pairs']))

if __name__=='__main__': unittest.main()
