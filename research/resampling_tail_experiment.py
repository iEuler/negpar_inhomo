"""Audit core/tail retention at fixed weight and with initial weight coarsening."""
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

NAMES=['mass','momentum_x','momentum_y','momentum_z','v2','anisotropy','fourth_x','fourth_radial']
for mode in ('100','010','001','110','200','020'):
    NAMES.extend((f'fourier_{mode}_real',f'fourier_{mode}_imag'))
CUTOFFS=(0.,2.5,3.,4.)
LABELS={0.:'Full wrapped',2.5:'Core 2.5 sigma',3.:'Core 3 sigma',4.:'Core 4 sigma'}

def errors(rows):
    return np.column_stack([rows[f'sample_{j}']-rows[f'source_{j}'] for j in range(len(NAMES))])

def analyze(data,replicas,rounds):
    expected=2*3*4*replicas*rounds
    if len(data)!=expected or not all(np.isfinite(data[n]).all() for n in data.dtype.names):
        raise RuntimeError('Incomplete or nonfinite probe output')
    exact=(data['weight_ratio']==1)|(data['round']>1)
    if not np.all(data['tail_moment_error'][exact]==0):
        raise RuntimeError('Equal-weight tail was not retained exactly')
    for j in range(len(NAMES)):
        if not np.all(data[f'source_{j}']==data[f'source_{j}'][0]):
            raise RuntimeError('Frozen source differs between rows')
    result={'rows':[],'paired':[]}
    rng=np.random.default_rng(831004)
    boot=rng.integers(0,replicas,size=(2000,replicas))
    for ratio in (1,4):
        for frequency in (4,8,12):
            for round_ in range(1,rounds+1):
                selected={}
                for cutoff in CUTOFFS:
                    d=data[(data['weight_ratio']==ratio)&(data['frequency']==frequency)&(data['round']==round_)&(data['cutoff']==cutoff)]
                    d=np.sort(d,order='replica')
                    if len(d)!=replicas or not np.array_equal(d['replica'],np.arange(replicas)):
                        raise RuntimeError('Missing or duplicate replica')
                    selected[cutoff]=d
                    delta=errors(d)
                    row=dict(weight_ratio=ratio,frequency=frequency,round=round_,cutoff=cutoff,
                        signed_count=float(np.mean(d['positive']+d['negative'])),
                        tail_count=float(np.mean(d['tail_positive']+d['tail_negative'])),
                        mean_seconds=float(d['seconds'].mean()),mean_cumulative_seconds=float(d['cumulative_seconds'].mean()),gate_accept_fraction=float(d['gate_accept'].mean()),
                        mean_attempts=float(d['attempts'].mean()),moments={})
                    for j,name in enumerate(NAMES):
                        row['moments'][name]=dict(mean_change=float(delta[:,j].mean()),
                            se=float(delta[:,j].std(ddof=1)/np.sqrt(replicas)),
                            rmse=float(np.sqrt(np.mean(delta[:,j]**2))))
                    row['fourier_rmse']=float(np.sqrt(np.mean(np.sum(delta[:,8:]**2,axis=1))))
                    result['rows'].append(row)
                if round_ in {1,rounds}:
                    base=selected[0.]; baseError=errors(base)
                    for cutoff in CUTOFFS[1:]:
                        d=selected[cutoff]; delta=errors(d)
                        pair=dict(weight_ratio=ratio,frequency=frequency,round=round_,cutoff=cutoff,rmse_ratios={},moment_changes={})
                        for name,j in (('anisotropy',5),('fourth_radial',7),('fourier',None)):
                            a=delta[:,j]**2 if j is not None else np.sum(delta[:,8:]**2,axis=1)
                            b=baseError[:,j]**2 if j is not None else np.sum(baseError[:,8:]**2,axis=1)
                            ratios=np.sqrt(a[boot].mean(axis=1)/b[boot].mean(axis=1))
                            pair['rmse_ratios'][name]=dict(ratio=float(np.sqrt(a.mean()/b.mean())),
                                ci95=[float(x) for x in np.quantile(ratios,[.025,.975])])
                        for j,name in enumerate(NAMES):
                            change=delta[:,j]-baseError[:,j]
                            pair['moment_changes'][name]=dict(mean=float(change.mean()),se=float(change.std(ddof=1)/np.sqrt(replicas)))
                        pair['count_ratio']=float(np.mean(d['positive']+d['negative'])/np.mean(base['positive']+base['negative']))
                        pair['time_ratio']=float(d['seconds'].mean()/base['seconds'].mean())
                        result['paired'].append(pair)
    return result

def write_report(summary,path):
    lines=['# Wrapped resampling with core/tail retention','',
        f"{summary['replicas']} replicas per mode/cutoff, {summary['rounds']} consecutive resampling calls. The source is the same frozen Gaussian-difference population as the geometry audit: 512 particles per sign, epsilon=0.05. All modes use certified quadratic envelopes and periodic wrapping. Fourier cutoffs are 4, 8, 12; physical core radii are 2.5, 3, 4 Maxwellian thermal units. The fixed background is zero mean, unit temperature.",'',
        'Weight ratio 1 keeps the original particle weight and retains tail particles exactly. Ratio 4 increases every particle weight fourfold on the first call; tails are uniformly thinned using the production partial-resampling rule, then retained exactly on subsequent calls. Counts include both signs and all retained tails. No moment correction, collision evolution, position reassignment, or weighted Fourier coupling is used. Timings include partition, reconstruction, thinning and merge; audit-only observable calculations are excluded. Mode order rotates and modes share seed schedules, but their proposal sequences/grids diverge, so pairing does not isolate identical proposals.','',
        'Every candidate is applied, even if it would fail the production requirement that both sign counts decrease. This intentionally exposes accumulated reconstruction error. Gate fractions are hypothetical; these are not solver trajectories with rejection/rollback. All errors are relative to the original frozen empirical source, not an exact evolving solution. Changing the core also changes its normalization box and effective physical resolution, so improvements cannot be attributed solely to the retained tail.','',
        '## Initial source tails','',
        '| Radius | Tail count | Share of unsigned radial fourth moment |','|---:|---:|---:|']
    for r in summary['source_tails']:
        lines.append(f"| {r['cutoff']} | {r['count']} | {100*r['unsigned_fourth_share']:.1f}% |")
    for ratio in (1,4):
        lines.extend(['',f'## Weight ratio {ratio}: cutoff 8','',
            '| Round | Method | Signed count | Call time (ms) | Cumulative time (ms) | Gate accepts | Mass drift | v² drift | Anisotropy RMSE | Radial fourth RMSE | Fourier RMSE |',
            '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|'])
        for r in summary['rows']:
            if r['weight_ratio']==ratio and r['frequency']==8 and r['round'] in {1,summary['rounds']}:
                m=r['moments']
                lines.append(f"| {r['round']} | {LABELS[r['cutoff']]} | {r['signed_count']:.1f} | {1000*r['mean_seconds']:.2f} | {1000*r['mean_cumulative_seconds']:.2f} | {100*r['gate_accept_fraction']:.0f}% | {m['mass']['mean_change']:.4g} | {m['v2']['mean_change']:.4g} | {m['anisotropy']['rmse']:.4g} | {m['fourth_radial']['rmse']:.4g} | {r['fourier_rmse']:.4g} |")
    lines.extend(['','## Paired error ratios against full wrapped resampling','',
        'Ratios below 1 indicate smaller RMSE. Parentheses show 95% paired bootstrap intervals (2,000 resamples). Extra retained particles can improve accuracy at additional cost; these are not matched-count or matched-accuracy efficiency claims.','',
        '| Weight ratio | Fourier cutoff | Round | Core radius | Count ratio | Time ratio | Anisotropy RMSE ratio (95% CI) | Radial fourth ratio (95% CI) | Fourier ratio (95% CI) |',
        '|---:|---:|---:|---:|---:|---:|---|---|---|'])
    for p in summary['paired']:
        fmt=lambda name: f"{p['rmse_ratios'][name]['ratio']:.3f} ({p['rmse_ratios'][name]['ci95'][0]:.3f}, {p['rmse_ratios'][name]['ci95'][1]:.3f})"
        lines.append(f"| {p['weight_ratio']} | {p['frequency']} | {p['round']} | {p['cutoff']} | {p['count_ratio']:.3f} | {p['time_ratio']:.3f} | {fmt('anisotropy')} | {fmt('fourth_radial')} | {fmt('fourier')} |")
    lines.extend(['','All 20 observables, standard errors, per-round metrics and paired comparisons are in summary.json. Raw samples, initial source velocities, console log, source snapshot and executable hash are preserved. Production defaults were not changed.'])
    (path/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def plot(summary,path):
    fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for i,ratio in enumerate((1,4)):
        for cutoff in CUTOFFS:
            r=[x for x in summary['rows'] if x['weight_ratio']==ratio and x['frequency']==8 and x['cutoff']==cutoff]
            rounds=[x['round'] for x in r]
            for j,values in enumerate(([x['signed_count'] for x in r],
                [x['moments']['anisotropy']['rmse'] for x in r],
                [x['moments']['fourth_radial']['rmse'] for x in r])):
                axes[i,j].plot(rounds,values,'o-',label=LABELS[cutoff])
        for j,label in enumerate(('Total signed particles','Anisotropy RMSE','Radial fourth-moment RMSE')):
            ax=axes[i,j]; ax.set_title(f'Weight ratio {ratio}: {label}'); ax.set_xlabel('Resampling call'); ax.grid(alpha=.2)
            ax.set_xticks(range(1,summary['rounds']+1)); ax.legend(fontsize=8)
    fig.suptitle('Repeated wrapped resampling at Fourier cutoff 8 — error relative to frozen input')
    fig.savefig(path/'repeated.png',dpi=180); plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for i,ratio in enumerate((1,4)):
        for cutoff in CUTOFFS:
            r=[x for x in summary['rows'] if x['weight_ratio']==ratio and x['round']==1 and x['cutoff']==cutoff]
            count=[x['signed_count'] for x in r]
            for j,values in enumerate(([x['moments']['anisotropy']['rmse'] for x in r],
                [x['moments']['fourth_radial']['rmse'] for x in r],[x['fourier_rmse'] for x in r])):
                axes[i,j].plot(count,values,'o-',label=LABELS[cutoff])
                for x,y,row in zip(count,values,r): axes[i,j].annotate(str(row['frequency']),(x,y),xytext=(4,4),textcoords='offset points',fontsize=7)
        for j,label in enumerate(('Anisotropy RMSE','Radial fourth-moment RMSE','Low Fourier-mode RMSE')):
            ax=axes[i,j]; ax.set_title(f'Weight ratio {ratio}: {label}'); ax.set_xlabel('Total signed particles after first call'); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('One-call accuracy versus retained count — labels mark Fourier cutoff')
    fig.savefig(path/'tradeoffs.png',dpi=180); plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replicas',type=int,default=128)
    parser.add_argument('--rounds',type=int,default=5)
    args=parser.parse_args()
    if args.replicas<8 or args.rounds<1: parser.error('Need replicas>=8 and rounds>=1')
    args.output.mkdir(parents=True)
    provenance=args.output/'provenance'; provenance.mkdir()
    executable=provenance/args.executable.name
    shutil.copy2(args.executable,executable)
    for dll in args.executable.parent.glob('*.dll'): shutil.copy2(dll,provenance/dll.name)
    root=Path(__file__).resolve().parents[1]
    with zipfile.ZipFile(provenance/'source.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for folder in ('src','research','tests'):
            for p in (root/folder).rglob('*'):
                if p.is_file() and 'runs' not in p.parts and '__pycache__' not in p.parts:
                    archive.write(p,p.relative_to(root))
        for name in ('CMakeLists.txt','negpar_inhomo.vcxproj'): archive.write(root/name,name)
    summary=dict(status='running',replicas=args.replicas,rounds=args.rounds,
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
    summary_path=args.output/'summary.json'
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    with (args.output/'console.log').open('w') as log:
        subprocess.run([str(executable.resolve()),str((args.output/'samples.csv').resolve()),str(args.replicas),str(args.rounds)],
            check=True,stdout=log,stderr=subprocess.STDOUT)
    data=np.genfromtxt(args.output/'samples.csv',delimiter=',',names=True)
    summary.update(analyze(data,args.replicas,args.rounds))
    initial=np.genfromtxt(args.output/'source_particles.csv',delimiter=',',names=True)
    radius2=initial['vx']**2+initial['vy']**2+initial['vz']**2
    summary['source_tails']=[dict(cutoff=c,count=int(np.sum(radius2>c*c)),
        unsigned_fourth_share=float(np.sum(radius2[radius2>c*c]**2)/np.sum(radius2**2))) for c in CUTOFFS[1:]]
    write_report(summary,args.output)
    plot(summary,args.output)
    summary['status']='complete'
    summary['samples_sha256']=hashlib.sha256((args.output/'samples.csv').read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    print(args.output/'REPORT.md')

if __name__=='__main__': main()
