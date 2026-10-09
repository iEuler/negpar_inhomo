"""Paired bounded signed-moment correction audit with aligned stratified cores."""
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


def analyze_corrected(data,replicas,rounds):
    if 'corrected' not in data.dtype.names or not np.all(np.isin(data['corrected'],[0,1])):
        raise RuntimeError('Invalid corrected method flags')
    if not np.all(data['aligned'][data['cutoff']>0]==1) or not np.all(data['corrected'][data['cutoff']==0]==0):
        raise RuntimeError('Study requires aligned core pairs and one independent full control')
    if not np.all(data['stratified'][data['cutoff']>0]==1):
        raise RuntimeError('Study requires stratified core pairs')
    successes=data[(data['corrected']==1)&(data['correction_status']==0)]
    for j in range(7):
        tolerance=1e-10*(successes['positive']+successes['negative'])*(.05/512)*successes['weight_ratio']*np.maximum(1.,successes['cutoff']**2)+1e-12
        if not np.all(np.abs(successes[f'correction_moment_{j}']-successes[f'call_target_{j}'])<=tolerance):
            raise RuntimeError('Successful correction violates measured moment constraint')
    mapped=data.copy(); mapped['aligned']=mapped['corrected']
    result=analyze_domains(mapped,replicas,rounds)
    result=dict(uncorrected=result['extrema'],corrected=result['aligned'],correction_pairs=result['geometry_pairs'])
    for method,flag in (('uncorrected',0),('corrected',1)):
        for row in result[method]['rows']:
            d=data[(data['weight_ratio']==row['weight_ratio'])&(data['frequency']==row['frequency'])&(data['cutoff']==row['cutoff'])&(data['round']==row['round'])]
            if row['cutoff']>0: d=d[d['corrected']==flag]
            row['attempts_sd']=float(d['attempts'].std(ddof=1))
            row['count_sd']=float((d['positive']+d['negative']).std(ddof=1))
            row['correction_success_fraction']=float(np.mean(d['correction_status']==0))
            row['correction_status_counts']={str(int(k)):int(v) for k,v in zip(*np.unique(d['correction_status'],return_counts=True))}
            successes=d[d['correction_status']==0]
            row['mean_applied_rms_displacement']=float(np.mean(successes['correction_rms'])) if len(successes) else 0.
            row['max_applied_displacement']=float(np.max(successes['correction_max'])) if len(successes) else 0.
            row['mean_attempted_removals']=float(np.mean(d['correction_removed']))
            row['mean_iterations']=float(np.mean(d['correction_iterations']))
            row['max_success_residual']=float(np.max(successes['correction_residual'])) if len(successes) else 0.
    return result


def verify_controls(data,previous):
    old=np.genfromtxt(previous,delimiter=',',names=True)
    a=data[data['corrected']==0]
    b=old[(old['stratified']==1)|(old['cutoff']==0)]
    keys=['weight_ratio','frequency','cutoff','replica','round']
    a=np.sort(a,order=keys); b=np.sort(b,order=keys)
    if len(a)!=len(b): raise RuntimeError('Control archive size mismatch')
    names=[n for n in b.dtype.names if n not in ('seconds','cumulative_seconds')]
    for n in names:
        if not np.array_equal(a[n],b[n]): raise RuntimeError('Previous control mismatch: '+n)
    return len(a)


def report(s,path):
    lines=['# Bounded signed-moment correction','',
        f"{s['replicas']} replicas, {s['rounds']} forced calls, Fourier cutoffs 4/8/12, core radii 2.5/3/4, output-weight ratios 1/4. Seven methods: full independent control, three aligned stratified cores, three aligned stratified cores with correction. Frozen Gaussian-difference source: 512 particles/sign, epsilon=0.05, seed 610701. Shared call seeds and rotating method order.",'',
        'Correction targets the input signed mass, momentum and three diagonal second moments of the entire group. Tail coarsening happens first; retained tail particles are fixed, and their output moments are subtracted from the total target to obtain the core target. Source mass is computed from signed counts. Equal-weight mass must be representable; excess-sign particles are uniformly deleted. The velocity solver takes minimum-norm linearized moment steps, freezes support-blocking particles, and backtracks within the physical sphere. It is a local iterative correction, not a global displacement optimum.','',
        'Limits fixed before this study: 40 iterations, normalized moment residual tolerance 1e-10 times retained count, RMS displacement at most 0.1 times core radius. Failure leaves the uncorrected core candidate unchanged. The audit then applies that uncorrected fallback at the requested output weight; it does not pretend an unsuccessful correction conserved moments. This differs from production rejection/rollback of the whole resampling candidate. Removal/displacement diagnostics on failures describe attempted changes only; applied displacement statistics use successful corrections.','',
        'Every resampling candidate is applied, including those failing the production count gate. No collisions or position reassignment. Equal-weight tails are exact; 4x weights thin tails on the first call. Errors are relative to the original frozen empirical source, not a continuum solution. Unconditional RMSE includes correction failures. Constrained moment accuracy alone is not evidence of distribution improvement: fourth moments and nonzero Fourier modes remain unconstrained. Timing includes target extraction, partition, reconstruction, tail adjustment, correction and merge, but excludes audit-only observables.','',
        f"All {s['control_reproduction_rows']} uncorrected control rows reproduce the preceding stratification archive in every non-timing field. Source snapshot, executable/DLLs, raw data, initial particles and hashes are archived. Production defaults are unchanged.",'']
    for ratio in (1,4):
        lines += [f'## Weight ratio {ratio}: radius 3, Fourier cutoff 8','',
            '| Call | Method | Count | Cumulative ms | Gate accepts | Correction succeeds | Applied RMS displacement | Mass RMSE | v2 RMSE | Anisotropy RMSE | Radial fourth RMSE | Fourier RMSE |',
            '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for round_ in sorted({1,s['rounds']}):
            for method in ('uncorrected','corrected'):
                r=next(x for x in s[method]['rows'] if x['weight_ratio']==ratio and x['frequency']==8 and x['cutoff']==3 and x['round']==round_)
                m=r['moments']
                lines.append(f"| {round_} | {method} | {r['signed_count']:.1f} | {1000*r['mean_cumulative_seconds']:.2f} | {100*r['gate_accept_fraction']:.0f}% | {100*r['correction_success_fraction']:.0f}% | {r['mean_applied_rms_displacement']:.4g} | {m['mass']['rmse']:.5g} | {m['v2']['rmse']:.5g} | {m['anisotropy']['rmse']:.5g} | {m['fourth_radial']['rmse']:.5g} | {r['fourier_rmse']:.5g} |")
        lines.append('')
    lines += ['## Corrected / uncorrected RMSE ratios','',
        'Below one favors correction. Parentheses: pointwise 95% paired bootstrap intervals, 2,000 resamples, conditional on this one source. Pairing shares seeds but source populations diverge after the first call. Mean/sign-count changes and correction feasibility are part of the method.','',
        '| Weight ratio | Fourier cutoff | Call | Radius | Count ratio | Cumulative time ratio | Mass | v2 | Anisotropy | Radial fourth | Fourier |',
        '|---:|---:|---:|---:|---:|---:|---|---|---|---|---|']
    for r in s['correction_pairs']:
        fmt=lambda n:f"{r['rmse_ratios'][n]['ratio']:.3f} ({r['rmse_ratios'][n]['ci95'][0]:.3f}, {r['rmse_ratios'][n]['ci95'][1]:.3f})"
        lines.append(f"| {r['weight_ratio']} | {r['frequency']} | {r['round']} | {r['cutoff']} | {r['count_ratio']:.3f} | {r['cumulative_time_ratio']:.3f} | {fmt('mass')} | {fmt('v2')} | {fmt('anisotropy')} | {fmt('fourth_radial')} | {fmt('fourier')} |")
    lines += ['','All 20 observables, mean changes and SEs, all per-call correction statuses, displacement/iteration/count diagnostics, and comparisons with the shared full control are retained in summary.json.']
    (path/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def plots(s,path):
    fig,axes=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
    for i,ratio in enumerate((1,4)):
        for method in ('uncorrected','corrected'):
            rows=[r for r in s[method]['rows'] if r['weight_ratio']==ratio and r['frequency']==8 and r['cutoff']==3]
            calls=[r['round'] for r in rows]
            for j,name in enumerate(('count','anisotropy','fourth_radial','fourier')):
                values=[r['signed_count'] if name=='count' else r['fourier_rmse'] if name=='fourier' else r['moments'][name]['rmse'] for r in rows]
                axes[i,j].plot(calls,values,'o-',label=method)
        for j,title in enumerate(('Total signed particles','Anisotropy RMSE','Radial fourth RMSE','Low Fourier RMSE')):
            axes[i,j].set_title(f'Weight ratio {ratio}: {title}',fontsize=10)
            axes[i,j].set_xlabel('Resampling call'); axes[i,j].set_xticks(range(1,s['rounds']+1)); axes[i,j].grid(alpha=.2)
            if j==0: axes[i,j].legend(fontsize=8)
    fig.suptitle('Aligned core radius 3, Fourier cutoff 8: bounded moment correction')
    fig.savefig(path/'repeated.png',dpi=180); plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
    for i,ratio in enumerate((1,4)):
        for round_ in sorted({1,s['rounds']}):
            rows=[r for r in s['correction_pairs'] if r['weight_ratio']==ratio and r['frequency']==8 and r['round']==round_]
            for j,name in enumerate(('anisotropy','fourth_radial','fourier')):
                value=np.array([r['rmse_ratios'][name]['ratio'] for r in rows]); ci=np.array([r['rmse_ratios'][name]['ci95'] for r in rows]).T
                axes[i,j].errorbar([r['cutoff'] for r in rows],value,yerr=np.maximum(0,np.vstack([value-ci[0],ci[1]-value])),fmt='o-',capsize=3,label=f'Call {round_}')
        for j,title in enumerate(('Anisotropy','Radial fourth','Low Fourier modes')):
            ax=axes[i,j]; ax.set_title(f'Weight ratio {ratio}: {title}'); ax.set_xlabel('Core radius (thermal units)'); ax.set_ylabel('Corrected / uncorrected RMSE'); ax.set_xticks([2.5,3,4]); ax.axhline(1,color='black',ls='--',lw=.8); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Bounded moment correction at Fourier cutoff 8: paired 95% bootstrap intervals')
    fig.savefig(path/'rmse_ratios.png',dpi=180); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replicas',type=int,default=128)
    parser.add_argument('--rounds',type=int,default=5)
    parser.add_argument('--previous',type=Path,default=Path(__file__).parent/'runs/resampling_stratified_v1/samples.csv')
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
        subprocess.run([str(executable.resolve()),str((args.output/'samples.csv').resolve()),str(args.replicas),str(args.rounds),'corrected'],check=True,stdout=log,stderr=subprocess.STDOUT)
    data=np.genfromtxt(args.output/'samples.csv',delimiter=',',names=True)
    s.update(analyze_corrected(data,args.replicas,args.rounds))
    s['control_reproduction_rows']=verify_controls(data,args.previous)
    report(s,args.output); plots(s,args.output)
    s['status']='complete'; s['samples_sha256']=hashlib.sha256((args.output/'samples.csv').read_bytes()).hexdigest()
    s['source_particles_sha256']=hashlib.sha256((args.output/'source_particles.csv').read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(s,indent=2)+'\n'); print(args.output/'REPORT.md')

if __name__=='__main__': main()
