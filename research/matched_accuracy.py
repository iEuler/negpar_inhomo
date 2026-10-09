"""Measured homogeneous cost at matched error, including complete HDP sampling.

No resampling. Tests ordinary PIC and count-proxy HDP mixing, plus a temporal
control variate using known initial moments for both methods. Selection is
over measured configurations, with independently validated selected points.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import zipfile
import numpy as np

OBSERVABLES=('anisotropy','fourth_moment','cosine_difference')


def initial(epsilon):
    return np.array([1.5*epsilon,3+9*epsilon,epsilon*(np.exp(-1)-np.exp(-.25))])


def scales(epsilon):
    return np.array([1.5*epsilon,9*epsilon,epsilon*(np.exp(-.25)-np.exp(-1))])


def read_run(directory):
    settings=json.loads((directory/'input.json').read_text())
    data=np.atleast_1d(np.genfromtxt(directory/'trajectories.csv',delimiter=',',names=True))
    count=settings['replicas'];steps=settings['steps']
    if len(data)!=count*(steps+1) or not all(np.isfinite(data[n]).all() for n in data.dtype.names):
        raise RuntimeError(f"Incomplete or nonfinite output: {directory}")
    cube=lambda prefix:np.stack([data[f'{prefix}_{j}'].reshape(count,steps+1) for j in range(3)],axis=-1)
    estimates={p:cube(p) for p in ('full','signed','mixed','pic_cv','mixed_cv')}
    if not np.all(data['full_count']==settings['full_count']):raise RuntimeError('Full count changed')
    for n in ('positive_count','negative_count'):
        if np.any(np.diff(data[n].reshape(count,steps+1),axis=1)<0):raise RuntimeError('Signed population decreased')
    costs=sum(data[n].reshape(count,steps+1)[:,-1] for n in
              ('initialization_seconds','collision_seconds','diagnostic_seconds'))
    return dict(estimates=estimates,costs=costs,settings=settings,
        count_max=int(np.max(data['positive_count']+data['negative_count'])),directory=str(directory))


def run(executable,directory,settings):
    if directory.exists():
        recorded=json.loads((directory/'input.json').read_text())
        if recorded!=settings:raise RuntimeError(f"Resume settings differ: {directory}")
        return read_run(directory)
    directory.mkdir()
    (directory/'input.json').write_text(json.dumps(settings,indent=2)+'\n')
    command=[str(executable),'--output',str((directory/'trajectories.csv').resolve())]
    for key,value in settings.items():command += ['--'+key.replace('_','-'),str(value)]
    start=time.perf_counter()
    with (directory/'console.log').open('w') as log:
        result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
    wall=time.perf_counter()-start
    (directory/'wall_seconds.json').write_text(json.dumps(wall)+'\n')
    if result.returncode:raise RuntimeError(f"Run failed; see {directory/'console.log'}")
    return read_run(directory)


def reference_estimates(run_data,epsilon):
    f=run_data['estimates']['full']
    # The known initial expectation reduces reference sampling uncertainty.
    # It changes the reference estimator, not the ordinary PIC competitor.
    return initial(epsilon)+f-f[:,:1,:]


def metrics(data,reference,epsilon,seed,bootstraps=300):
    truth=reference.mean(axis=0)
    normalization=scales(epsilon)
    rng=np.random.default_rng(seed)
    methods=('full','pic_cv') if data['settings']['mode']=='pic' else ('signed','mixed','mixed_cv')
    result={}
    bootstrap_indices=[rng.integers(len(data['costs']),size=len(data['costs'])) for _ in range(bootstraps)]
    ref_means=[reference[rng.integers(len(reference),size=len(reference))].mean(axis=0) for _ in range(bootstraps)]
    for method in methods:
        values=data['estimates'][method]
        errors=(values-truth)/normalization
        losses=(errors*errors).mean(axis=1) # per trajectory, per observable
        joint=losses.mean(axis=1)
        bootstrap=[]
        for indices,ref in zip(bootstrap_indices,ref_means):
            bootstrap.append(np.sqrt(np.mean(((values[indices]-ref)/normalization)**2)))
        result[method]=dict(joint_relative_rmse=float(np.sqrt(joint.mean())),
            joint_relative_rmse_ci95=[float(v) for v in np.quantile(bootstrap,[.025,.975])],
            observable_relative_rmse=[float(v) for v in np.sqrt(losses.mean(axis=0))],
            final_reference_bias=[float(v) for v in (values[:,-1,:].mean(axis=0)-truth[-1])],
            mean_compute_seconds=float(data['costs'].mean()),compute_seconds_se=float(data['costs'].std(ddof=1)/np.sqrt(len(data['costs']))))
    return result


def choose(rows,method,target):
    feasible=[r for r in rows if method in r['metrics'] and r['metrics'][method]['joint_relative_rmse_ci95'][1]<=target]
    return min(feasible,key=lambda r:r['metrics'][method]['mean_compute_seconds']) if feasible else None


def report(summary,output):
    lines=['# Complete homogeneous HDP mixing versus PIC at matched accuracy','',
        'Collision-only evolution; no resampling, synchronization, projection or advection. Timings include initialization, source-bound construction, full and signed collision/source evolution, all observable estimates, mixing, and temporal-control-variate arithmetic. Research CSV serialization is excluded from compute time; wall time including output is retained for each run.','',
        'Error is trajectory RMS over the three observables, each normalized by its known initial departure from Maxwellian equilibrium. Reference ensembles use their known initial expectation as a control variate; particle-count and half-timestep reference checks quantify sensitivity.','',
        'PIC and HDP are tuned over measured count configurations. Candidate selection requires the upper 95% bootstrap error bound to meet the target. Selected configurations are rerun with fresh seeds for validation. This is a measured search, not a proof of globally optimal configurations.','',
        'Ordinary comparison: PIC full estimate versus HDP count-proxy mixture. Stronger comparison: both methods use known initial moments as a temporal control variate. These do not require fitted pilot weights.','',
        '| epsilon | Target relative RMS | Comparison | Validated PIC RMS | Validated HDP RMS | PIC/HDP compute time (s) | PIC time / HDP time |','|---:|---:|---|---:|---:|---:|---:|']
    for case in summary['cases']:
        for pair in case.get('comparisons',[]):
            if not pair.get('complete'):
                lines.append(f"| {case['epsilon']:g} | {pair['target']:g} | {pair['comparison']} | unavailable | unavailable | unavailable | unavailable |")
                continue
            p,h=pair['pic'],pair['hdp']
            label=pair['comparison'] if pair['target_verified'] else pair['comparison']+' (target not verified)'
            lines.append(f"| {case['epsilon']:g} | {pair['target']:g} | {label} | {p['joint_relative_rmse']:.3f} | {h['joint_relative_rmse']:.3f} | {p['mean_compute_seconds']:.4g}/{h['mean_compute_seconds']:.4g} | {pair['time_ratio']:.3f} |")
    lines += ['', 'A ratio above one favors HDP, subject to both independently validated errors meeting the target. Different achieved errors are reported explicitly; this is a common-accuracy-threshold comparison rather than exact equality of errors.','', '## Reference checks','']
    for case in summary['cases']:
        lines.append(f"epsilon={case['epsilon']:g}: normalized reference mean uncertainty {case['reference_uncertainty']:.4g}; higher-count reference trajectory difference {case['reference_count_difference']:.4g}; half-dt reference difference {case['reference_dt_difference']:.4g}.")
    lines += ['', 'All candidate measurements, bootstrap error bounds, selected particle counts, fresh validation runs and provenance are saved in summary.json and the per-run directories. Bias is included in the measured error. The approximate kernel, finite radial support and empirical remainder envelope remain numerical limitations; this benchmark does not certify consistency or long-time conservation of the discretization.']
    if (output/'confirmation'/'REPORT.md').exists():
        lines += ['', '## Supplemental fresh confirmation', '',
                  'See [the closer-error confirmation](confirmation/REPORT.md) for additional independent runs, including a passing larger signed ensemble after the original epsilon=0.01 CV validation missed its bound. Original measurements and failed validation remain above.']
    (output/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--epsilons',type=float,nargs='+',default=[.05,.01,.3])
    parser.add_argument('--replicas',type=int,default=192)
    parser.add_argument('--large-replicas',type=int,default=64)
    parser.add_argument('--validation-replicas',type=int,default=384)
    parser.add_argument('--reference-replicas',type=int,default=128)
    parser.add_argument('--reference-count',type=int,default=32768)
    parser.add_argument('--sign-counts',type=int,nargs='+',default=[16,32,64,128,256])
    parser.add_argument('--full-ratios',type=int,nargs='+',default=[8,16])
    parser.add_argument('--pic-counts',type=int,nargs='+',default=[32,128,512,2048,8192,32768,131072,524288])
    parser.add_argument('--targets',type=float,nargs='+',default=[.25,.4])
    parser.add_argument('--steps',type=int,default=20)
    parser.add_argument('--dt',type=float,default=.01)
    parser.add_argument('--strength',type=float,default=5.)
    args=parser.parse_args()
    if min(args.replicas,args.large_replicas,args.validation_replicas,args.reference_replicas)<8:parser.error('Use at least eight replicas per ensemble')
    args.output.mkdir(parents=True,exist_ok=args.resume)
    provenance=args.output/'provenance'
    executable=provenance/args.executable.name
    if not args.resume:
        provenance.mkdir()
        shutil.copy2(args.executable,executable)
        for dll in args.executable.parent.glob('*.dll'):shutil.copy2(dll,provenance/dll.name)
        root=Path(__file__).resolve().parents[1]
        with zipfile.ZipFile(provenance/'source.zip','w',zipfile.ZIP_DEFLATED) as archive:
            for folder in ('src','research'):
                for path in (root/folder).rglob('*'):
                    if path.is_file() and 'runs' not in path.parts and '__pycache__' not in path.parts:archive.write(path,path.relative_to(root))
            archive.write(root/'CMakeLists.txt','CMakeLists.txt')
    executable=executable.resolve()
    summary=dict(status='running',arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),cases=[])
    summary_path=args.output/'summary.json'
    try:
        for ci,epsilon in enumerate(args.epsilons):
            directory=args.output/f'epsilon_{epsilon:g}'
            directory.mkdir(exist_ok=args.resume)
            base=700000000+ci*300000000
            def settings(mode,nf,np_,replicas,seed,dt=args.dt,steps=args.steps):
                return dict(mode=mode,full_count=nf,sign_count=np_,replicas=replicas,seed=seed,epsilon=epsilon,dt=dt,steps=steps,strength=args.strength)
            print(f'epsilon={epsilon} references',flush=True)
            r=run(executable,directory/'reference',settings('pic',args.reference_count,16,args.reference_replicas,base))
            reference=reference_estimates(r,epsilon)
            check_count=max(32,args.reference_replicas//4)
            large=run(executable,directory/'reference_count_double',settings('pic',2*args.reference_count,16,check_count,base+10000000))
            refined=run(executable,directory/'reference_dt_half',settings('pic',args.reference_count,16,check_count,base+20000000,args.dt/2,2*args.steps))
            norm=scales(epsilon)
            truth=reference.mean(axis=0)
            case=dict(epsilon=epsilon,candidates=[],comparisons=[],
                reference_uncertainty=float(np.sqrt(np.mean((reference.std(axis=0,ddof=1)/np.sqrt(len(reference))/norm)**2))),
                reference_count_difference=float(np.sqrt(np.mean(((reference_estimates(large,epsilon).mean(axis=0)-truth)/norm)**2))),
                reference_dt_difference=float(np.sqrt(np.mean(((reference_estimates(refined,epsilon)[:,::2,:].mean(axis=0)-truth)/norm)**2))))
            summary['cases'].append(case)
            def candidate(label,config):
                print(f'epsilon={epsilon} {label}',flush=True)
                measured=run(executable,directory/label,config)
                row=dict(label=label,settings=config,metrics=metrics(measured,reference,epsilon,config['seed']+99),
                    signed_count_max=measured['count_max'],directory=str(directory/label))
                case['candidates'].append(row)
                summary_path.write_text(json.dumps(summary,indent=2)+'\n')
                return row
            for ni,np_ in enumerate(args.sign_counts):
                for ri,ratio in enumerate(args.full_ratios):
                    candidate(f'hdp_n{np_}_ratio{ratio}',settings('hdp',np_*ratio,np_,args.replicas,base+30000000+ni*2000000+ri*1000000))
            for ni,nf in enumerate(args.pic_counts):
                count=args.replicas if nf<=32768 else args.large_replicas
                candidate(f'pic_n{nf}',settings('pic',nf,16,count,base+60000000+ni*1000000))
                if choose(case['candidates'],'full',min(args.targets)) and choose(case['candidates'],'pic_cv',min(args.targets)):break
            # Refine count brackets by the empirical affine MSE versus 1/N
            # model; measure every proposed point rather than claim extrapolation.
            for method in ('full','pic_cv'):
                for ti,target in enumerate(args.targets):
                    rows=[x for x in case['candidates'] if method in x['metrics']]
                    below=[x for x in rows if x['metrics'][method]['joint_relative_rmse']>target]
                    above=[x for x in rows if x['metrics'][method]['joint_relative_rmse']<=target]
                    if not below or not above:continue
                    low=max(below,key=lambda x:x['settings']['full_count'])
                    high=min(above,key=lambda x:x['settings']['full_count'])
                    n1,n2=low['settings']['full_count'],high['settings']['full_count']
                    if n1>=n2:continue
                    e1,e2=(x['metrics'][method]['joint_relative_rmse']**2 for x in (low,high))
                    slope=(e1-e2)/(1/n1-1/n2);intercept=e2-slope/n2
                    if slope<=0 or target*target<=intercept:continue
                    estimate=int(np.ceil(1.2*slope/(target*target-intercept)/32))*32
                    if n1<estimate<n2 and all(abs(estimate-x['settings']['full_count'])>.08*estimate for x in rows):
                        candidate(f'pic_refine_{method}_target{target:g}_n{estimate}',settings('pic',estimate,16,args.replicas if estimate<=32768 else args.large_replicas,base+80000000+ti*1000000+(5000000 if method=='pic_cv' else 0)))
            validation_cache={}
            for ti,target in enumerate(args.targets):
                for pi,(p_method,h_method,name) in enumerate((('full','mixed','ordinary'),('pic_cv','mixed_cv','initial-moment CV'))):
                    selected=[choose(case['candidates'],m,target) for m in (p_method,h_method)]
                    pair=dict(target=target,comparison=name,complete=False)
                    case['comparisons'].append(pair)
                    if any(x is None for x in selected):continue
                    validated=[]
                    for si,(chosen,method) in enumerate(zip(selected,(p_method,h_method))):
                        key=chosen['label']
                        if key not in validation_cache:
                            config=dict(chosen['settings'])
                            config.update(seed=base+120000000+len(validation_cache)*2000000,
                                replicas=args.validation_replicas if config['full_count']<=32768 else max(64,args.large_replicas))
                            measured=run(executable,directory/f'validate_{key}',config)
                            validation_cache[key]=dict(settings=config,metrics=metrics(measured,reference,epsilon,config['seed']+77))
                        validated.append(validation_cache[key]['metrics'][method])
                    p,h=validated
                    pair.update(complete=True,pic=p,hdp=h,pic_settings=validation_cache[selected[0]['label']]['settings'],
                        hdp_settings=validation_cache[selected[1]['label']]['settings'],
                        target_verified=all(x['joint_relative_rmse_ci95'][1]<=target for x in validated),
                        time_ratio=p['mean_compute_seconds']/h['mean_compute_seconds'])
                    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
            report(summary,args.output)
        summary['status']='complete'
    except Exception as error:
        summary.update(status='failed',failure=str(error));raise
    finally:summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    report(summary,args.output)
    print(args.output/'REPORT.md')


if __name__=='__main__':main()
