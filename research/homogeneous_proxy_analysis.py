"""Bootstrap final-time gains of the existing population-count mixing proxy."""
import argparse
import json
from pathlib import Path
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('experiment',type=Path)
    args=parser.parse_args()
    summary=json.loads((args.experiment/'summary.json').read_text())
    rows=[]
    for index,case in enumerate(summary['cases']):
        folder=args.experiment/f"epsilon_{case['epsilon']:g}"
        data=np.genfromtxt(folder/'evaluation/trajectories.csv',delimiter=',',names=True)
        ref=np.genfromtxt(folder/'reference/trajectories.csv',delimiter=',',names=True)
        final=data[data['step']==summary['arguments']['steps']]
        reference=ref[ref['step']==summary['arguments']['steps']]
        f=np.column_stack([final[f'full_{j}'] for j in range(3)])
        d=np.column_stack([final[f'signed_{j}'] for j in range(3)])
        r=np.column_stack([reference[f'full_{j}'] for j in range(3)])
        w=final['count_proxy_weight'][:,None]
        mixed=w*f+(1-w)*d
        rng=np.random.default_rng(71000+index)
        ratios=[]
        for _ in range(800):
            selection=rng.integers(len(f),size=len(f))
            truth=r[rng.integers(len(r),size=len(r))].mean(axis=0)
            mf=np.mean((f[selection]-truth)**2,axis=0)
            md=np.mean((d[selection]-truth)**2,axis=0)
            mm=np.mean((mixed[selection]-truth)**2,axis=0)
            ratios.append(np.minimum(mf,md)/mm)
        ci=np.quantile(ratios,[.025,.975],axis=0)
        truth=r.mean(axis=0)
        gain=np.minimum(np.mean((f-truth)**2,axis=0),np.mean((d-truth)**2,axis=0))/np.mean((mixed-truth)**2,axis=0)
        for j,name in enumerate(('anisotropy','fourth_moment','cosine_difference')):
            rows.append(dict(epsilon=case['epsilon'],observable=name,mse_gain=float(gain[j]),
                             ci95=[float(v) for v in ci[:,j]],mean_full_weight=float(w.mean())))
    lines=['# Existing count-proxy mixing: homogeneous final-time test','',
        'This uses the same evaluation trajectories and independent references as homogeneous_v1. The weight uses only current particle counts and weights, with no pilot calibration. Bootstrap intervals jointly resample evaluation and reference replicas. Gain > 1 favors mixing. These are individual intervals, not simultaneous confidence statements.','',
        '| epsilon | Observable | MSE gain (95% interval) | Mean full weight |','|---:|---|---:|---:|']
    for row in rows:
        lo,hi=row['ci95']
        lines.append(f"| {row['epsilon']:g} | {row['observable']} | {row['mse_gain']:.3f} [{lo:.3f}, {hi:.3f}] | {row['mean_full_weight']:.3f} |")
    (args.experiment/'PROXY_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (args.experiment/'proxy_summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    (args.experiment/'provenance/proxy_analysis_source.py').write_bytes(Path(__file__).read_bytes())
    print(args.experiment/'PROXY_REPORT.md')


if __name__=='__main__':main()
