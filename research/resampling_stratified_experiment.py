"""Paired stratified proposal-allocation audit with aligned physical cores."""
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
from resampling_domain_experiment import analyze_domains, NAMES


def analyze_stratified(data,replicas,rounds):
    if 'stratified' not in data.dtype.names or not np.all(np.isin(data['stratified'],[0,1])):
        raise RuntimeError('Invalid stratified method flags')
    if not np.all(data['aligned'][data['cutoff']>0]==1) or not np.all(data['stratified'][data['cutoff']==0]==0):
        raise RuntimeError('Study requires aligned core pairs and one independent full control')
    mapped=data.copy(); mapped['aligned']=mapped['stratified']
    result=analyze_domains(mapped,replicas,rounds)
    result=dict(independent=result['extrema'],stratified=result['aligned'],allocation_pairs=result['geometry_pairs'])
    for method,flag in (('independent',0),('stratified',1)):
        for row in result[method]['rows']:
            d=data[(data['weight_ratio']==row['weight_ratio'])&(data['frequency']==row['frequency'])&(data['cutoff']==row['cutoff'])&(data['round']==row['round'])]
            if row['cutoff']>0: d=d[d['stratified']==flag]
            row['attempts_sd']=float(d['attempts'].std(ddof=1))
            row['count_sd']=float((d['positive']+d['negative']).std(ddof=1))
    return result


def verify_controls(data,previous):
    old=np.genfromtxt(previous,delimiter=',',names=True)
    a=data[data['stratified']==0]
    b=old[(old['aligned']==1)|(old['cutoff']==0)]
    keys=['weight_ratio','frequency','cutoff','replica','round']
    a=np.sort(a,order=keys); b=np.sort(b,order=keys)
    if len(a)!=len(b): raise RuntimeError('Control archive size mismatch')
    names=[n for n in b.dtype.names if n not in ('seconds','cumulative_seconds')]
    for n in names:
        if not np.array_equal(a[n],b[n]): raise RuntimeError('Previous control mismatch: '+n)
    return len(a)


def report(s,path):
    lines=['# Stratified proposal allocation with aligned physical cores','',
        f"{s['replicas']} replicas, {s['rounds']} forced calls, Fourier cutoffs 4/8/12, radii 2.5/3/4, output-weight ratios 1/4. Seven methods: full independent control, three aligned independent cores, three aligned stratified cores. Frozen Gaussian-difference source: 512 particles/sign, epsilon=0.05, seed 610701. Calls share seed schedules; method order rotates.",'',
        'Stratification partitions cumulative expected envelope proposals into unit intervals, with an independent uniform point in each interval. Cells follow lexicographic grid order. Only boundary strata need explicit draws; fully covered strata contribute deterministically. Every cell retains its expected proposal count. The total is floor/ceiling of the cumulative expectation. Uniform positions and independent rejection decisions are unchanged. This is proposal-allocation stratification, not spatial stratification within a cell and not a general variance guarantee for signed observables. Certified envelopes are required; production defaults are unchanged.','',
        'Both core arms use the same fixed physical bounds, partition, weights, Fourier/Taylor reconstruction, sphere mask and tail adjustment. Equal-weight tails are retained exactly; 4x weights thin tails on the first call. All candidates are applied, including those failing the production reduction gate. No moment correction, collisions, or position reassignment. Errors are relative to one frozen empirical source, not a continuum solution. Timing includes partition/reconstruction/thinning/merge, excludes audit observables; extra particles imply future solver cost not measured here.','',
        f"All {s['control_reproduction_rows']} independent control rows reproduce the preceding domain archive exactly for every non-timing field. Raw data, source particles, source snapshot, executable/DLLs and hashes are archived.",'']
    for ratio in (1,4):
        lines += [f'## Weight ratio {ratio}: radius 3, Fourier cutoff 8','',
            '| Call | Allocation | Count | Count SD | Attempts SD | Cumulative ms | Gate accepts | Mass drift +/- SE | v2 drift +/- SE | Anisotropy RMSE | Radial fourth RMSE | Fourier RMSE |',
            '|---:|---|---:|---:|---:|---:|---:|---|---|---:|---:|---:|']
        for round_ in sorted({1,s['rounds']}):
            for method in ('independent','stratified'):
                r=next(x for x in s[method]['rows'] if x['weight_ratio']==ratio and x['frequency']==8 and x['cutoff']==3 and x['round']==round_)
                m=r['moments']; fmt=lambda n:f"{m[n]['mean_change']:.5g} +/- {m[n]['se']:.3g}"
                lines.append(f"| {round_} | {method} | {r['signed_count']:.1f} | {r['count_sd']:.2f} | {r['attempts_sd']:.2f} | {1000*r['mean_cumulative_seconds']:.2f} | {100*r['gate_accept_fraction']:.0f}% | {fmt('mass')} | {fmt('v2')} | {m['anisotropy']['rmse']:.5g} | {m['fourth_radial']['rmse']:.5g} | {r['fourier_rmse']:.5g} |")
        lines.append('')
    lines += ['## Stratified / independent RMSE ratios','',
        'Below one favors stratification. Parentheses: pointwise 95% paired bootstrap intervals, 2,000 resamples, conditional on this frozen source; not simultaneous sweep intervals. Shared seeds pair replicas but changed RNG consumption means proposal locations differ.','',
        '| Weight ratio | Fourier cutoff | Call | Radius | Count ratio | Cumulative time ratio | Mass | v2 | Anisotropy | Radial fourth | Fourier |',
        '|---:|---:|---:|---:|---:|---:|---|---|---|---|---|']
    for r in s['allocation_pairs']:
        fmt=lambda n:f"{r['rmse_ratios'][n]['ratio']:.3f} ({r['rmse_ratios'][n]['ci95'][0]:.3f}, {r['rmse_ratios'][n]['ci95'][1]:.3f})"
        lines.append(f"| {r['weight_ratio']} | {r['frequency']} | {r['round']} | {r['cutoff']} | {r['count_ratio']:.3f} | {r['cumulative_time_ratio']:.3f} | {fmt('mass')} | {fmt('v2')} | {fmt('anisotropy')} | {fmt('fourth_radial')} | {fmt('fourier')} |")
    lines += ['','summary.json includes all 20 observables for every call and setting, mean errors and standard errors, counts, proposal/count SD, runtimes, gate rates, and paired RMSE ratios.']
    (path/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def plots(s,path):
    fig,axes=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
    for i,ratio in enumerate((1,4)):
        for method in ('independent','stratified'):
            rows=[r for r in s[method]['rows'] if r['weight_ratio']==ratio and r['frequency']==8 and r['cutoff']==3]
            calls=[r['round'] for r in rows]
            for j,name in enumerate(('count','anisotropy','fourth_radial','fourier')):
                values=[r['signed_count'] if name=='count' else r['fourier_rmse'] if name=='fourier' else r['moments'][name]['rmse'] for r in rows]
                axes[i,j].plot(calls,values,'o-',label=method)
        for j,title in enumerate(('Total signed particles','Anisotropy RMSE','Radial fourth RMSE','Low Fourier RMSE')):
            axes[i,j].set_title(f'Weight ratio {ratio}: {title}',fontsize=10)
            axes[i,j].set_xlabel('Resampling call'); axes[i,j].set_xticks(range(1,s['rounds']+1)); axes[i,j].grid(alpha=.2)
            if j==0: axes[i,j].legend(fontsize=8)
    fig.suptitle('Aligned core radius 3, Fourier cutoff 8: proposal-allocation comparison')
    fig.savefig(path/'repeated.png',dpi=180); plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
    for i,ratio in enumerate((1,4)):
        for round_ in sorted({1,s['rounds']}):
            rows=[r for r in s['allocation_pairs'] if r['weight_ratio']==ratio and r['frequency']==8 and r['round']==round_]
            for j,name in enumerate(('anisotropy','fourth_radial','fourier')):
                value=np.array([r['rmse_ratios'][name]['ratio'] for r in rows]); ci=np.array([r['rmse_ratios'][name]['ci95'] for r in rows]).T
                axes[i,j].errorbar([r['cutoff'] for r in rows],value,yerr=np.maximum(0,np.vstack([value-ci[0],ci[1]-value])),fmt='o-',capsize=3,label=f'Call {round_}')
        for j,title in enumerate(('Anisotropy','Radial fourth','Low Fourier modes')):
            ax=axes[i,j]; ax.set_title(f'Weight ratio {ratio}: {title}'); ax.set_xlabel('Core radius (thermal units)'); ax.set_ylabel('Stratified / independent RMSE'); ax.set_xticks([2.5,3,4]); ax.axhline(1,color='black',ls='--',lw=.8); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Proposal stratification at Fourier cutoff 8: paired 95% bootstrap intervals')
    fig.savefig(path/'rmse_ratios.png',dpi=180); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replicas',type=int,default=128)
    parser.add_argument('--rounds',type=int,default=5)
    parser.add_argument('--previous',type=Path,default=Path(__file__).parent/'runs/resampling_domain_v1/samples.csv')
    args=parser.parse_args()
    if args.replicas<8 or not 1<=args.rounds<=100: parser.error('Need replicas>=8 and rounds in 1..100')
    args.output.mkdir(parents=True); provenance=args.output/'provenance'; provenance.mkdir()
    executable=provenance/args.executable.name; shutil.copy2(args.executable,executable)
    for dll in args.executable.parent.glob('*.dll'): shutil.copy2(dll,provenance/dll.name)
    root=Path(__file__).resolve().parents[1]
    with zipfile.ZipFile(provenance/'source.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for folder in ('src','research','tests'):
            for p in (root/folder).rglob('*'):
                if p.is_file() and 'runs' not in p.parts and '__pycache__' not in p.parts: archive.write(p,p.relative_to(root))
        for name in ('CMakeLists.txt','negpar_inhomo.vcxproj'): archive.write(root/name,name)
    s=dict(status='running',replicas=args.replicas,rounds=args.rounds,executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
    summary_path=args.output/'summary.json'; summary_path.write_text(json.dumps(s,indent=2)+'\n')
    with (args.output/'console.log').open('w') as log:
        subprocess.run([str(executable.resolve()),str((args.output/'samples.csv').resolve()),str(args.replicas),str(args.rounds),'stratified'],check=True,stdout=log,stderr=subprocess.STDOUT)
    data=np.genfromtxt(args.output/'samples.csv',delimiter=',',names=True)
    s.update(analyze_stratified(data,args.replicas,args.rounds))
    s['control_reproduction_rows']=verify_controls(data,args.previous)
    report(s,args.output); plots(s,args.output)
    s['status']='complete'; s['samples_sha256']=hashlib.sha256((args.output/'samples.csv').read_bytes()).hexdigest()
    s['source_particles_sha256']=hashlib.sha256((args.output/'source_particles.csv').read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(s,indent=2)+'\n'); print(args.output/'REPORT.md')

if __name__=='__main__': main()
