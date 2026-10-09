import unittest
import numpy as np
from test_resampling_tail_experiment import fixture
from resampling_domain_experiment import analyze_domains, NAMES

def domain_fixture():
    old=fixture()
    aligned=old[old['cutoff']>0].copy()
    for j in range(len(NAMES)):
        aligned[f'sample_{j}']=aligned[f'source_{j}']+.5*(aligned[f'sample_{j}']-aligned[f'source_{j}'])
    data=np.zeros(len(old)+len(aligned),dtype=old.dtype.descr+[('aligned','f8')])
    for name in old.dtype.names:
        data[name][:len(old)]=old[name]; data[name][len(old):]=aligned[name]
    data['aligned'][len(old):]=1
    return data

class DomainTests(unittest.TestCase):
    def test_geometry_pairs_and_shared_control(self):
        result=analyze_domains(domain_fixture()[::-1],8,2)
        self.assertEqual(len(result['extrema']['rows']),48)
        self.assertEqual(len(result['aligned']['rows']),48)
        for pair in result['geometry_pairs']:
            for value in pair['rmse_ratios'].values():
                self.assertAlmostEqual(value['ratio'],.5)
                self.assertAlmostEqual(value['ci95'][0],.5)
                self.assertAlmostEqual(value['ci95'][1],.5)

    def test_geometry_source_mismatch_rejected(self):
        data=domain_fixture(); data['source_0'][-1]=1
        with self.assertRaisesRegex(RuntimeError,'Source mismatch'):
            analyze_domains(data,8,2)

if __name__=='__main__': unittest.main()
