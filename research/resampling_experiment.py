"""Isolate the resampling envelope change on one frozen signed population."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile
import numpy as np

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replicas',type=int,default=512)
    parser.add_argument('--geometry',action='store_true',help='Include certified periodic wrapping')
    args=parser.parse_args()
    args.output.mkdir(parents=True)
    provenance=args.output/'provenance'; provenance.mkdir()
    executable=provenance/args.executable.name
    shutil.copy2(args.executable,executable)
    for dll in args.executable.parent.glob('*.dll'): shutil.copy2(dll,provenance/dll.name)
    root=Path(__file__).resolve().parents[1]
    with zipfile.ZipFile(provenance/'source.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for folder in ('src','research','tests'):
            for path in (root/folder).rglob('*'):
                if path.is_file() and 'runs' not in path.parts and '__pycache__' not in path.parts:
                    archive.write(path,path.relative_to(root))
        archive.write(root/'CMakeLists.txt','CMakeLists.txt')
    summary=dict(status='running',replicas=args.replicas,geometry=args.geometry,
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),rows=[],paired_changes=[])
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (args.output/'console.log').open('w') as log:
        subprocess.run([str(executable.resolve()),str((args.output/'samples.csv').resolve()),str(args.replicas)]+(['geometry'] if args.geometry else []),
                       check=True,stdout=log,stderr=subprocess.STDOUT)
    data=np.genfromtxt(args.output/'samples.csv',delimiter=',',names=True)
    modes=[(0,0,'Legacy adaptive'),(1,0,'Certified shifted')]+([(1,1,'Certified wrapped')] if args.geometry else [])
    if len(data)!=args.replicas*3*len(modes) or not all(np.isfinite(data[n]).all() for n in data.dtype.names):
        raise RuntimeError('Incomplete or nonfinite resampling results')
    names=('mass','momentum_x','v2','anisotropy','fourth_moment','cosine_difference')
    for frequency in (4,8,12):
        for certified,wrapped,label in modes:
            d=data[(data['frequency']==frequency)&(data['certified']==certified)&(data['wrapped']==wrapped)]
            row=dict(frequency=frequency,certified=bool(certified),wrapped=bool(wrapped),label=label,
                mean_seconds=float(d['seconds'].mean()),mean_attempts=float(d['attempts'].mean()),
                total_envelope_increases=int(d['envelope_increases'].sum()),
                mean_signed_count=float((d['positive']+d['negative']).mean()),moments={})
            for j,name in enumerate(names):
                delta=d[f'sample_{j}']-d[f'source_{j}']
                row['moments'][name]=dict(mean_change=float(delta.mean()),
                    se=float(delta.std(ddof=1)/np.sqrt(len(delta))),rmse=float(np.sqrt(np.mean(delta**2))))
            summary['rows'].append(row)
        old=data[(data['frequency']==frequency)&(data['certified']==0)]
        new=data[(data['frequency']==frequency)&(data['certified']==1)&(data['wrapped']==0)]
        paired={}
        for j,name in enumerate(names):
            delta=new[f'sample_{j}']-old[f'sample_{j}']
            paired[name]=dict(mean=float(delta.mean()),se=float(delta.std(ddof=1)/np.sqrt(len(delta))))
        summary['paired_changes'].append(dict(frequency=frequency,moments=paired))
        if args.geometry:
            wrapped=data[(data['frequency']==frequency)&(data['wrapped']==1)]
            assert np.array_equal(new['replica'],wrapped['replica'])
            assert np.array_equal(new['attempts'],wrapped['attempts'])
            changes={}
            for j,name in enumerate(names):
                delta=wrapped[f'sample_{j}']-new[f'sample_{j}']
                changes[name]=dict(mean=float(delta.mean()),se=float(delta.std(ddof=1)/np.sqrt(len(delta))))
            summary.setdefault('geometry_paired_changes',[]).append(dict(frequency=frequency,moments=changes))
    summary['status']='complete'
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Fixed-envelope Fourier resampling audit','',
        f"{args.replicas} replicas per mode and cutoff, one frozen signed Gaussian-difference population (512 of each sign), epsilon=0.05. Output weight is four times input weight. Coefficients, quadratic interpolation and Fourier cutoffs are shared. The certified wrapped mode periodically maps proposal coordinates into [0,2*pi) before sphere masking and velocity restoration; its Taylor offsets and fixed envelope remain those of the original local cell. No moment correction or collision evolution is used. Timing order rotates within each replica and includes the full reconstruction.",'',
        'The certified bound is for the piecewise quadratic reconstruction, not the original distribution or exact Fourier interpolant. Remaining truncation, interpolation, spherical support and sampling noise affect moments. Source-to-output RMSE is a diagnostic relative to the frozen particle estimator, not physical accuracy against an evolving solution.','',
        '| Cutoff | Envelope | Time (ms) | Attempts | Bound increases | Signed count | Mass change ± SE | v² change ± SE | Anisotropy RMSE |',
        '|---:|---|---:|---:|---:|---:|---|---|---:|']
    for r in summary['rows']:
        m=r['moments']; fmt=lambda n:f"{m[n]['mean_change']:.5g} ± {m[n]['se']:.3g}"
        lines.append(f"| {r['frequency']} | {r['label']} | {1000*r['mean_seconds']:.3f} | {r['mean_attempts']:.1f} | {r['total_envelope_increases']} | {r['mean_signed_count']:.1f} | {fmt('mass')} | {fmt('v2')} | {m['anisotropy']['rmse']:.5g} |")
    if args.geometry:
        lines+=['','## Paired wrapping changes (wrapped minus shifted certified)', '',
            'Both certified modes use identical proposals and acceptance draws; the differences below isolate added spherical-cap samples. Errors relative to the input still include Fourier/Taylor and sphere-truncation effects.', '',
            '| Cutoff | Mass change ± SE | v² change ± SE |','|---:|---|---|']
        for r in summary['geometry_paired_changes']:
            fmt=lambda n:f"{r['moments'][n]['mean']:.6g} ± {r['moments'][n]['se']:.3g}"
            lines.append(f"| {r['frequency']} | {fmt('mass')} | {fmt('v2')} |")
    lines+=['','All six moment changes, RMSEs and paired differences are retained in summary.json. This does not validate repeated resampling in Landau damping or adaptive synchronization rollback. The default solver remains on the legacy mode.']
    (args.output/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for certified,wrapped,label in modes:
        rows=[r for r in summary['rows'] if r['certified']==bool(certified) and r['wrapped']==bool(wrapped)]
        axes[0].plot([r['frequency'] for r in rows],[r['mean_seconds']*1000 for r in rows],'o-',label=label)
        axes[1].plot([r['frequency'] for r in rows],[r['moments']['anisotropy']['rmse'] for r in rows],'o-',label=label)
    for ax in axes: ax.set_xlabel('Fourier cutoff'); ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[0].set_ylabel('Reconstruction time (ms)');axes[1].set_ylabel('Anisotropy RMSE vs frozen input')
    fig.savefig(args.output/'comparison.png',dpi=180)
    print(args.output/'REPORT.md')

if __name__=='__main__': main()
