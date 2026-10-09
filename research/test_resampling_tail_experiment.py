"""Check paired statistics and rejection of incomplete/corrupted audit data."""
import unittest
import numpy as np
from resampling_tail_experiment import analyze, NAMES, CUTOFFS

def fixture():
    fields=['weight_ratio','frequency','cutoff','replica','round','seconds','cumulative_seconds',
        'attempts','positive','negative','tail_positive','tail_negative','gate_accept','tail_moment_error']
    for j in range(len(NAMES)): fields.extend((f'source_{j}',f'sample_{j}'))
    data=np.zeros(2*3*4*8*2,dtype=[(name,'f8') for name in fields])
    index=0
    for ratio in (1,4):
        for frequency in (4,8,12):
            for round_ in (1,2):
                for cutoff in CUTOFFS:
                    for replica in range(8):
                        row=data[index]; index+=1
                        row['weight_ratio']=ratio; row['frequency']=frequency
                        row['cutoff']=cutoff; row['replica']=replica; row['round']=round_
                        row['positive']=20; row['negative']=20
                        row['seconds']=.01; row['cumulative_seconds']=round_*.01
                        scale=1. if cutoff==0 else .5
                        for j in range(len(NAMES)):
                            row[f'source_{j}']=j*.1
                            row[f'sample_{j}']=j*.1+scale*(replica+1)*.001
    return data

class AuditTests(unittest.TestCase):
    def test_paired_error_ratio_and_bootstrap(self):
        data=fixture()[::-1]  # analysis must pair by replica, not incoming row order
        result=analyze(data,8,2)
        for pair in result['paired']:
            for value in pair['rmse_ratios'].values():
                self.assertAlmostEqual(value['ratio'],.5)
                self.assertAlmostEqual(value['ci95'][0],.5)
                self.assertAlmostEqual(value['ci95'][1],.5)

    def test_duplicate_replica_is_rejected(self):
        data=fixture(); data['replica'][1]=data['replica'][0]
        with self.assertRaisesRegex(RuntimeError,'Missing or duplicate'):
            analyze(data,8,2)

    def test_retention_error_only_allowed_during_weight_change(self):
        data=fixture()
        coarsened=(data['weight_ratio']==4)&(data['round']==1)
        data['tail_moment_error'][coarsened]=.1
        analyze(data,8,2)
        data['tail_moment_error'][0]=.1
        with self.assertRaisesRegex(RuntimeError,'retained exactly'):
            analyze(data,8,2)

if __name__=='__main__': unittest.main()
