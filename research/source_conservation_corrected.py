"""Fresh trajectories after fixing the source rejection count dependency."""
import json
from pathlib import Path
import shutil
import numpy as np
from matched_accuracy import run, read_run, metrics, reference_estimates


def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'research/runs/source_conservation_corrected_v1'
    out.mkdir(exist_ok=True)
    exe=out/'negpar_homogeneous.exe'
    if not exe.exists():
        shutil.copy2(root/'build/release/Release/negpar_homogeneous.exe',exe)
        for dll in (root/'build/release/Release').glob('*.dll'):shutil.copy2(dll,out/dll.name)
        shutil.copy2(root/'src/resampling/NegativeParticleSampling.cpp',out/'sampling.cpp')
        shutil.copy2(root/'research/homogeneous_mixing.cpp',out/'harness.cpp')
        shutil.copy2(Path(__file__),out/'analysis.py')
    summary=[]
    configs=[('eps0.3',.3,.01,20,2048,128,1024),
             ('eps0.3_halfdt',.3,.005,40,2048,128,1024),
             ('eps0.05',.05,.01,20,2048,128,1024),
             ('eps0.05_halfdt',.05,.005,40,2048,128,1024),
             ('eps0.01_cv',.01,.01,20,256,32,768),
             ('eps0.01_ordinary',.01,.01,20,1024,128,384)]
    for i,(label,epsilon,dt,steps,nf,np_,replicas) in enumerate(configs):
        config=dict(mode='hdp',full_count=nf,sign_count=np_,replicas=replicas,
                    seed=2500000000+i*1000000,epsilon=epsilon,dt=dt,steps=steps,strength=5.)
        print(label,flush=True)
        data=run(exe.resolve(),out/label,config)
        raw=np.genfromtxt(out/label/'trajectories.csv',delimiter=',',names=True)
        row=dict(label=label,settings=config,drift={})
        for prefix in ('signed','full'):
            for name,j in [('mass',3),('px',4),('py',5),('pz',6),('total_v2',7)]:
                values=raw[f'{prefix}_{j}'].reshape(replicas,steps+1)
                drift=values[:,-1]-values[:,0]
                row['drift'][f'{prefix}_{name}']=dict(mean=float(drift.mean()),se=float(drift.std(ddof=1)/np.sqrt(replicas)))
        if epsilon==.01:
            archive=root/'research/runs/matched_accuracy_v1'
            reference=reference_estimates(read_run(archive/'epsilon_0.01/reference'),epsilon)
            row['metrics']=metrics(data,reference,epsilon,config['seed']+99)
            pic_label='pic_cv_close' if label.endswith('cv') else 'pic_ordinary_close'
            pic_method='pic_cv' if label.endswith('cv') else 'full'
            hdp_method='mixed_cv' if label.endswith('cv') else 'mixed'
            pic=json.loads((archive/'confirmation/summary.json').read_text())[pic_label]['metrics'][pic_method]
            row['pic_comparison']=dict(pic=pic,hdp_method=hdp_method,
                time_ratio=pic['mean_compute_seconds']/row['metrics'][hdp_method]['mean_compute_seconds'])
        summary.append(row)
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(row,flush=True)
    print(out/'summary.json',flush=True)


if __name__=='__main__':main()
