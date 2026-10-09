"""Validate source reconstruction and sampling after the exact split fix."""
import json
from pathlib import Path
import shutil
import subprocess
import zipfile
import numpy as np
from matched_accuracy import run, read_run, metrics, reference_estimates


def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'research/runs/source_split_v1'
    out.mkdir(exist_ok=True)
    provenance=out/'provenance'
    if not provenance.exists():
        provenance.mkdir()
        for name in ('negpar_source_probe.exe','negpar_homogeneous.exe','libfftw3-3.dll'):
            shutil.copy2(root/'build/release/Release'/name,provenance/name)
        with zipfile.ZipFile(provenance/'source.zip','w',zipfile.ZIP_DEFLATED) as archive:
            for folder in ('src','research','tests'):
                for path in (root/folder).rglob('*'):
                    if path.is_file() and 'runs' not in path.parts and '__pycache__' not in path.parts:
                        archive.write(path,path.relative_to(root))
    single=out/'single_source.csv'
    if not single.exists():
        print('single-source sampling',flush=True)
        subprocess.run([str(provenance/'negpar_source_probe.exe'),'single',str(single),'.01','30000'],check=True)
    data=np.genfromtxt(single,delimiter=',',names=True)
    summary={'single_source':[],'trajectories':[]}
    for speed in (1.,2.):
        p=data[(data['speed']==speed)&(data['sign']==1)]
        n=data[(data['speed']==speed)&(data['sign']==-1)]
        # Swapping signed lists reverses their summation order. Mass is exact;
        # floating velocity moments can differ by final summation roundoff.
        error=max(float(np.max(np.abs(p[k]+n[k]))) for k in ('px','v2'))
        symmetry=bool(np.array_equal(p['mass'],-n['mass']) and error<1e-13)
        if not symmetry:raise RuntimeError('Signed source reversal failed')
        row=dict(speed=speed,replicas=len(p),mass_mean=float(p['mass'].mean()),
                 mass_se=float(p['mass'].std(ddof=1)/np.sqrt(len(p))),sign_reversal_to_roundoff=symmetry,
                 sign_reversal_max_moment_error=error)
        summary['single_source'].append(row);print(row,flush=True)
    for i,(epsilon,dt) in enumerate(((.3,.01),(.3,.005),(.05,.01),(.05,.005),(.01,.01))):
        label=f'epsilon{epsilon:g}_dt{dt:g}'
        print(label,flush=True)
        nf,np_,replicas=(1024,128,384) if epsilon==.01 else (2048,128,512)
        config=dict(mode='hdp',full_count=nf,sign_count=np_,replicas=replicas,
                    seed=3100000000+i*1000000,epsilon=epsilon,dt=dt,steps=round(.2/dt),strength=5.)
        measured=run((provenance/'negpar_homogeneous.exe').resolve(),out/label,config)
        raw=np.genfromtxt(out/label/'trajectories.csv',delimiter=',',names=True)
        row=dict(settings=config,drift={},compute_seconds=float(measured['costs'].mean()))
        if epsilon==.01:
            archive=root/'research/runs/matched_accuracy_v1'
            ref=reference_estimates(read_run(archive/'epsilon_0.01/reference'),epsilon)
            row['metrics']=metrics(measured,ref,epsilon,config['seed']+99)
            pic=json.loads((archive/'confirmation/summary.json').read_text())['pic_ordinary_close']['metrics']['full']
            row['pic_ordinary']=pic
            row['time_ratio']=pic['mean_compute_seconds']/row['compute_seconds']
        for name,j in (('mass',3),('px',4),('py',5),('pz',6),('v2',7)):
            values=raw[f'signed_{j}'].reshape(replicas,config['steps']+1)
            change=values[:,-1]-values[:,0]
            row['drift'][name]=dict(mean=float(change.mean()),se=float(change.std(ddof=1)/np.sqrt(replicas)))
        summary['trajectories'].append(row)
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(row,flush=True)
    lines=['# Exact signed-source split validation','',
           'The bounded branch uses clamp(q,-a,a); its signed remainder is q minus that bounded piece. The bounded proposal and its count use the sum of both signed populations. Both proposal branches use the collision-background density. Out-of-bound acceptance ratios raise an error instead of silently clipping.', '',
           '## Single-source mass and sign reversal', '',
           '| Source speed | Replicas | Mean mass change ± SE | Sign reversal within roundoff |',
           '|---:|---:|---:|---|']
    for row in summary['single_source']:
        lines.append(f"| {row['speed']:g} | {row['replicas']} | {row['mass_mean']:.6g} ± {row['mass_se']:.3g} | {row['sign_reversal_to_roundoff']} |")
    lines += ['', '## No-resampling trajectories, T=0.2, A=5', '',
              '| epsilon | dt | Replicas | Signed mass drift ± SE | Signed total-v² drift ± SE |',
              '|---:|---:|---:|---:|---:|']
    for row in summary['trajectories']:
        s=row['settings'];m=row['drift']['mass'];e=row['drift']['v2']
        lines.append(f"| {s['epsilon']} | {s['dt']} | {s['replicas']} | {m['mean']:.6g} ± {m['se']:.3g} | {e['mean']:.6g} ± {e['se']:.3g} |")
    near=summary['trajectories'][-1]
    if 'metrics' in near:
        h=near['metrics']['mixed'];p=near['pic_ordinary']
        lines += ['', '## Selected near-equilibrium efficiency point', '',
                  f"Ordinary PIC RMS {p['joint_relative_rmse']:.3f} versus repaired HDP mixture RMS {h['joint_relative_rmse']:.3f} (95% interval {h['joint_relative_rmse_ci95']}). PIC compute {p['mean_compute_seconds']:.4g} s versus HDP {near['compute_seconds']:.4g} s; time ratio {near['time_ratio']:.2f}. Reuses unaffected pure PIC runs; this is a selected-point recheck, not a new optimization."]
    lines += ['', 'All trajectories completed with explicit rejection-bound checks enabled. This restores the source decomposition; kernel approximation, finite radial support and global certification of the empirical remainder bound remain separate limitations. No invariant projection or resampling was used. Frozen executables, source archive, inputs and raw samples are retained.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(out/'REPORT.md',flush=True)


if __name__=='__main__':main()
