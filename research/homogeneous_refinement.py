"""Check signed conservation drift after halving dt, still without resampling."""
import argparse
import json
from pathlib import Path
import numpy as np
from homogeneous_experiment import run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--replicas",type=int,default=512)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    summary=json.loads((args.experiment/'summary.json').read_text())
    settings=summary['arguments']
    executable=(args.experiment/'provenance/negpar_homogeneous.exe').resolve()
    records=[]
    for index,case in enumerate(summary['cases']):
        epsilon=case['epsilon']
        print(f"epsilon={epsilon} signed dt/2",flush=True)
        refined=run(executable,args.output/f"epsilon_{epsilon:g}",replicas=args.replicas,
            steps=2*settings['steps'],full_count=settings['full_count'],sign_count=settings['sign_count'],
            epsilon=epsilon,dt=settings['dt']/2,strength=settings['strength'],seed=170000000+index*100000000)
        coarse=np.genfromtxt(args.experiment/f"epsilon_{epsilon:g}/evaluation/trajectories.csv",delimiter=',',names=True)
        for j,name in ((3,'mass'),(7,'total_second_moment')):
            values=coarse[f'signed_{j}'].reshape(settings['replicas'],settings['steps']+1)
            delta=values[:,-1]-values[:,0]
            fine_delta=refined['signed'][:,-1,j]-refined['signed'][:,0,j]
            records.append(dict(epsilon=epsilon,observable=name,coarse_drift=float(delta.mean()),
                coarse_drift_se=float(delta.std(ddof=1)/np.sqrt(len(delta))),
                half_dt_drift=float(fine_delta.mean()),half_dt_drift_se=float(fine_delta.std(ddof=1)/np.sqrt(len(fine_delta)))))
    lines=['# Homogeneous signed-conservation timestep sensitivity','',
        'Independent ensembles at the same final time; no resampling, projection or advection. Uncertainty is the standard error of each trajectory change, accounting for its initial/final correlation. This diagnoses the existing source sampler, not a convergence order.','',
        '| epsilon | Observable | Coarse drift ± SE | Half-dt drift ± SE |','|---:|---|---:|---:|']
    for r in records:
        lines.append(f"| {r['epsilon']:g} | {r['observable']} | {r['coarse_drift']:.5g} ± {r['coarse_drift_se']:.2g} | {r['half_dt_drift']:.5g} ± {r['half_dt_drift_se']:.2g} |")
    (args.output/'summary.json').write_text(json.dumps(records,indent=2)+'\n')
    (args.output/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (args.output/'analysis_source.py').write_bytes(Path(__file__).read_bytes())
    print(args.output/'REPORT.md')


if __name__=='__main__':main()
