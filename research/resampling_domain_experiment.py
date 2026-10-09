"""Paired audit of fixed physical core bounds versus per-axis particle extrema."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from resampling_tail_experiment import analyze, errors, NAMES, CUTOFFS, LABELS

def geometry_pairs(data,replicas,rounds):
    if len(data)!=2*3*7*replicas*rounds or not all(np.isfinite(data[n]).all() for n in data.dtype.names):
        raise RuntimeError('Incomplete or nonfinite domain audit')
    pairs=[]
    boot=np.random.default_rng(831005).integers(0,replicas,size=(2000,replicas))
    for ratio in (1,4):
        for frequency in (4,8,12):
            for round_ in sorted({1,rounds}):
                for cutoff in CUTOFFS[1:]:
                    selected=(data['weight_ratio']==ratio)&(data['frequency']==frequency)&(data['round']==round_)&(data['cutoff']==cutoff)
                    old=np.sort(data[selected&(data['aligned']==0)],order='replica')
                    new=np.sort(data[selected&(data['aligned']==1)],order='replica')
                    if len(old)!=replicas or len(new)!=replicas or not np.array_equal(old['replica'],np.arange(replicas)) or not np.array_equal(old['replica'],new['replica']):
                        raise RuntimeError('Missing or duplicate paired replica')
                    for j in range(len(NAMES)):
                        if not np.array_equal(old[f'source_{j}'],new[f'source_{j}']):
                            raise RuntimeError('Source mismatch between geometries')
                    a,b=errors(new),errors(old)
                    row=dict(weight_ratio=ratio,frequency=frequency,round=round_,cutoff=cutoff,rmse_ratios={},moment_changes={})
                    for name,j in (('mass',0),('v2',4),('anisotropy',5),('fourth_radial',7),('fourier',None)):
                        a2=a[:,j]**2 if j is not None else np.sum(a[:,8:]**2,axis=1)
                        b2=b[:,j]**2 if j is not None else np.sum(b[:,8:]**2,axis=1)
                        ratios=np.sqrt(a2[boot].mean(axis=1)/b2[boot].mean(axis=1))
                        row['rmse_ratios'][name]=dict(ratio=float(np.sqrt(a2.mean()/b2.mean())),
                            ci95=[float(x) for x in np.quantile(ratios,[.025,.975])])
                    for j,name in enumerate(NAMES):
                        d=a[:,j]-b[:,j]
                        row['moment_changes'][name]=dict(mean=float(d.mean()),se=float(d.std(ddof=1)/np.sqrt(replicas)))
                    row['count_ratio']=float(np.mean(new['positive']+new['negative'])/np.mean(old['positive']+old['negative']))
                    row['time_ratio']=float(new['seconds'].mean()/old['seconds'].mean())
                    row['cumulative_time_ratio']=float(new['cumulative_seconds'].mean()/old['cumulative_seconds'].mean())
                    pairs.append(row)
    return pairs

def analyze_domains(data,replicas,rounds):
    pairs=geometry_pairs(data,replicas,rounds)
    extrema=analyze(data[data['aligned']==0],replicas,rounds)
    aligned=analyze(data[(data['aligned']==1)|(data['cutoff']==0)],replicas,rounds)
    return dict(extrema=extrema,aligned=aligned,geometry_pairs=pairs)

def report(s,path):
    lines=['# Aligned physical core bounds','',
        f"{s['replicas']} replicas per method/setting and {s['rounds']} consecutive forced resampling calls. One frozen Gaussian-difference source (512 per sign, epsilon=0.05); Fourier cutoffs 4/8/12, core radii 2.5/3/4 thermal units, output-weight ratios 1 and 4. All seven methods run together with rotating order and shared seed schedules: full wrapped, three extrema-based core methods, and three fixed-bound core methods.",'',
        'Aligned core bounds are u_i +/- cutoff*sqrt(T), with zero mean and unit temperature held fixed here. The normalized spherical support then equals the physical core sphere. The control computes its bounds from each core population\'s axis extrema, so its support is a physical ellipsoid. The full wrapped control retains its extrema-based bounds. Both core methods use the same partition rule, certified quadratic envelopes, wrapping, Fourier reconstruction and tail adjustment. Grid spacing changes with the box, so this audit measures the complete domain change, not support masking alone.','',
        'Weight ratio 1 retains tails exactly; ratio 4 uniformly thins tails on the first call, then retains them exactly. No moment correction, collisions, or position reassignment. Every candidate is applied to expose reconstruction error accumulation; the production count gate is measured but not enforced. Counts include all tails. Timings include partition/reconstruction/thinning/merge and exclude audit-only observable calculations. All errors are relative to the original empirical source, not a continuum solution. Shared seeds pair replicas, but different grids/populations change proposal sequences.','']
    for ratio in (1,4):
        lines.extend([f'## Weight ratio {ratio}, Fourier cutoff 8, core radius 3','',
            '| Call | Method | Count | Cumulative ms | Gate accepts | Mass drift +/- SE | v² drift +/- SE | Anisotropy RMSE | Radial fourth RMSE | Fourier RMSE |',
            '|---:|---|---:|---:|---:|---|---|---:|---:|---:|'])
        for round_ in sorted({1,s['rounds']}):
            for geometry,cutoff,label in (('extrema',0,'Full wrapped'),('extrema',3,'Core extrema'),('aligned',3,'Core aligned')):
                r=next(x for x in s[geometry]['rows'] if x['weight_ratio']==ratio and x['frequency']==8 and x['cutoff']==cutoff and x['round']==round_)
                m=r['moments']; fmt=lambda n:f"{m[n]['mean_change']:.5g} +/- {m[n]['se']:.3g}"
                lines.append(f"| {round_} | {label} | {r['signed_count']:.1f} | {1000*r['mean_cumulative_seconds']:.2f} | {100*r['gate_accept_fraction']:.0f}% | {fmt('mass')} | {fmt('v2')} | {m['anisotropy']['rmse']:.5g} | {m['fourth_radial']['rmse']:.5g} | {r['fourier_rmse']:.5g} |")
        lines.append('')
    lines.extend(['## Aligned/extrema paired RMSE ratios','',
        'Below 1 favors aligned bounds; parentheses are pointwise 95% paired bootstrap intervals (2,000 resamples). They are conditional on this single source and are not simultaneous claims over the sweep.','',
        '| Weight ratio | Fourier cutoff | Call | Radius | Count ratio | Call-time ratio | Mass ratio (95% CI) | v² ratio (95% CI) | Anisotropy ratio (95% CI) | Radial fourth ratio (95% CI) | Fourier ratio (95% CI) |',
        '|---:|---:|---:|---:|---:|---:|---|---|---|---|---|'])
    for r in s['geometry_pairs']:
        fmt=lambda n:f"{r['rmse_ratios'][n]['ratio']:.3f} ({r['rmse_ratios'][n]['ci95'][0]:.3f}, {r['rmse_ratios'][n]['ci95'][1]:.3f})"
        lines.append(f"| {r['weight_ratio']} | {r['frequency']} | {r['round']} | {r['cutoff']} | {r['count_ratio']:.3f} | {r['time_ratio']:.3f} | {fmt('mass')} | {fmt('v2')} | {fmt('anisotropy')} | {fmt('fourth_radial')} | {fmt('fourier')} |")
    lines.extend(['','All 20 observable statistics, per-call count/cost/gate measurements, aligned/extrema pairs, and each method\'s comparisons with full wrapped resampling are in summary.json. Raw samples, initial particles, executable/DLLs, source ZIP and hashes are preserved. Production defaults are unchanged.'])
    (path/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def plots(s,path):
    fig,axes=plt.subplots(2,4,figsize=(16,7),constrained_layout=True)
    for i,ratio in enumerate((1,4)):
        for geometry,cutoff,label in (('extrema',0,'Full wrapped'),('extrema',3,'Core extrema 3 sigma'),('aligned',3,'Core aligned 3 sigma')):
            r=[x for x in s[geometry]['rows'] if x['weight_ratio']==ratio and x['frequency']==8 and x['cutoff']==cutoff]
            calls=[x['round'] for x in r]
            for j,name in enumerate(('mass','v2','anisotropy','fourier')):
                if j<2:
                    axes[i,j].errorbar(calls,[x['moments'][name]['mean_change'] for x in r],
                        yerr=[1.96*x['moments'][name]['se'] for x in r],fmt='o-',capsize=3,label=label)
                else:
                    values=[x['fourier_rmse'] if name=='fourier' else x['moments'][name]['rmse'] for x in r]
                    axes[i,j].plot(calls,values,'o-',label=label)
        for j,title in enumerate(('Mass drift (95% mean interval)','v² drift (95% mean interval)','Anisotropy RMSE','Low Fourier-mode RMSE')):
            ax=axes[i,j]; ax.set_title(f'Weight ratio {ratio}: {title}',fontsize=10); ax.set_xlabel('Resampling call'); ax.grid(alpha=.2)
            ax.set_xticks(range(1,s['rounds']+1))
            if j<2: ax.axhline(0,color='black',lw=.7,alpha=.5)
            if j==0: ax.legend(fontsize=7)
    fig.suptitle('Physical-domain comparison: radius 3, Fourier cutoff 8, no moment correction')
    fig.savefig(path/'domain_comparison.png',dpi=180); plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
    for i,ratio in enumerate((1,4)):
        for round_ in sorted({1,s['rounds']}):
            r=[x for x in s['geometry_pairs'] if x['weight_ratio']==ratio and x['frequency']==8 and x['round']==round_]
            for j,name in enumerate(('anisotropy','fourth_radial','fourier')):
                value=np.array([x['rmse_ratios'][name]['ratio'] for x in r])
                interval=np.array([x['rmse_ratios'][name]['ci95'] for x in r]).T
                axes[i,j].errorbar([x['cutoff'] for x in r],value,yerr=np.maximum(0,np.vstack([value-interval[0],interval[1]-value])),fmt='o-',capsize=3,label=f'Call {round_}')
        for j,title in enumerate(('Anisotropy','Radial fourth moment','Low Fourier modes')):
            ax=axes[i,j]; ax.set_title(f'Weight ratio {ratio}: {title}'); ax.set_xlabel('Core radius (thermal units)'); ax.set_ylabel('Aligned / extrema RMSE')
            ax.axhline(1,color='black',lw=.8,ls='--'); ax.grid(alpha=.2); ax.legend(fontsize=8); ax.set_xticks(CUTOFFS[1:])
    fig.suptitle('Aligned versus extrema bounds at Fourier cutoff 8 — paired 95% bootstrap intervals')
    fig.savefig(path/'rmse_ratios.png',dpi=180); plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replicas',type=int,default=128)
    parser.add_argument('--rounds',type=int,default=5)
    args=parser.parse_args()
    if args.replicas<8 or not 1<=args.rounds<=100: parser.error('Need replicas>=8 and rounds in 1..100')
    args.output.mkdir(parents=True)
    provenance=args.output/'provenance'; provenance.mkdir()
    executable=provenance/args.executable.name
    shutil.copy2(args.executable,executable)
    for dll in args.executable.parent.glob('*.dll'): shutil.copy2(dll,provenance/dll.name)
    root=Path(__file__).resolve().parents[1]
    with zipfile.ZipFile(provenance/'source.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for folder in ('src','research','tests'):
            for p in (root/folder).rglob('*'):
                if p.is_file() and 'runs' not in p.parts and '__pycache__' not in p.parts: archive.write(p,p.relative_to(root))
        for name in ('CMakeLists.txt','negpar_inhomo.vcxproj'): archive.write(root/name,name)
    summary=dict(status='running',replicas=args.replicas,rounds=args.rounds,
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
    summary_path=args.output/'summary.json'; summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    with (args.output/'console.log').open('w') as log:
        subprocess.run([str(executable.resolve()),str((args.output/'samples.csv').resolve()),str(args.replicas),str(args.rounds),'aligned'],check=True,stdout=log,stderr=subprocess.STDOUT)
    data=np.genfromtxt(args.output/'samples.csv',delimiter=',',names=True)
    summary.update(analyze_domains(data,args.replicas,args.rounds))
    report(summary,args.output); plots(summary,args.output)
    summary['status']='complete'; summary['samples_sha256']=hashlib.sha256((args.output/'samples.csv').read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    print(args.output/'REPORT.md')

if __name__=='__main__': main()
