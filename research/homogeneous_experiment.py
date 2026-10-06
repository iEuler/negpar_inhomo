"""Replicated homogeneous Coulomb mixing test, strictly without resampling."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import zipfile
import numpy as np

NAMES = ("anisotropy", "fourth_moment", "cosine_difference", "mass", "momentum_x", "momentum_y", "momentum_z", "total_second_moment")


def run(executable, directory, *, replicas, steps, full_count, sign_count, epsilon, dt, strength, seed, mode="hdp"):
    directory.mkdir()
    output=directory / "trajectories.csv"
    command=[str(executable), "--output", str(output.resolve())]
    settings=dict(replicas=replicas, steps=steps, full_count=full_count, sign_count=sign_count,
                  epsilon=epsilon, dt=dt, strength=strength, seed=seed, mode=mode)
    for key,value in settings.items():
        command += ["--"+key.replace("_","-"),str(value)]
    (directory / "input.json").write_text(json.dumps(settings,indent=2)+"\n")
    start=time.perf_counter()
    with (directory / "console.log").open("w") as log:
        result=subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
    if result.returncode:
        raise RuntimeError(f"Homogeneous run failed: {directory / 'console.log'}")
    wall=time.perf_counter()-start
    data=np.genfromtxt(output, delimiter=",", names=True)
    data=np.atleast_1d(data)
    if len(data)!=replicas*(steps+1) or not all(np.isfinite(data[name]).all() for name in data.dtype.names):
        raise RuntimeError(f"Invalid/incomplete homogeneous output: {output}")
    cube=lambda prefix: np.stack([data[f"{prefix}_{j}"].reshape(replicas,steps+1) for j in range(8)],axis=-1)
    full, signed=cube("full"),cube("signed")
    counts={name:data[name].reshape(replicas,steps+1) for name in ("positive_count","negative_count","full_count")}
    if not np.all(counts['full_count']==full_count):
        raise RuntimeError("Full population changed in a no-resampling run")
    for name in ("positive_count","negative_count"):
        if np.any(np.diff(counts[name],axis=1)<0):
            raise RuntimeError("Signed population decreased without a permitted removal operation")
    numerical={name:data[name].copy() for name in data.dtype.names if not name.endswith("seconds")}
    cost={name:float(np.mean(data[name].reshape(replicas,steps+1)[:,-1]))
          for name in ("initialization_seconds","collision_seconds","diagnostic_seconds")}
    return dict(full=full,signed=signed,counts=counts,proxy=data['count_proxy_weight'].reshape(replicas,steps+1),
                cost=cost,wall_seconds=wall,numerical=numerical)


def fit(pilot):
    f,d=pilot['full'],pilot['signed']
    vf,vd=f.var(axis=0,ddof=1),d.var(axis=0,ddof=1)
    covariance=((f-f.mean(axis=0))*(d-d.mean(axis=0))).sum(axis=0)/(len(f)-1)
    denominator=vf+vd-2*covariance
    weight=np.divide(vd-covariance,denominator,out=np.zeros_like(vd),where=denominator>1e-20)
    independent=np.divide(vd,vf+vd,out=np.zeros_like(vd),where=vf+vd>1e-20)
    return np.clip(weight,0,1),np.clip(independent,0,1),vf,vd,covariance


def bootstrap_gain(full,signed,mixed,reference,exact_initial,rng):
    ratios=[]
    for _ in range(400):
        indices=rng.integers(len(full),size=len(full))
        ref=reference[rng.integers(len(reference),size=len(reference))].mean(axis=0)
        ref[0]=exact_initial
        mf=np.mean((full[indices]-ref)**2,axis=0)
        md=np.mean((signed[indices]-ref)**2,axis=0)
        mm=np.mean((mixed[indices]-ref)**2,axis=0)
        ratios.append(np.minimum(mf,md)/np.maximum(mm,1e-30))
    return np.quantile(ratios,[.025,.975],axis=0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--pilot",type=int,default=256)
    parser.add_argument("--replicas",type=int,default=1024)
    parser.add_argument("--reference-replicas",type=int,default=64)
    parser.add_argument("--full-count",type=int,default=2048)
    parser.add_argument("--sign-count",type=int,default=128)
    parser.add_argument("--reference-count",type=int,default=32768)
    parser.add_argument("--steps",type=int,default=20)
    parser.add_argument("--dt",type=float,default=.01)
    parser.add_argument("--strength",type=float,default=5.)
    parser.add_argument("--epsilons",type=float,nargs="+",default=[.05,.3])
    args=parser.parse_args()
    if min(args.pilot,args.replicas,args.reference_replicas)<2:
        parser.error("Every ensemble needs at least two replicas")
    args.output.mkdir(parents=True,exist_ok=False)
    provenance=args.output / "provenance"
    provenance.mkdir()
    executable=provenance / args.executable.name
    shutil.copy2(args.executable,executable)
    for dll in args.executable.parent.glob("*.dll"):
        shutil.copy2(dll,provenance / dll.name)
    root=Path(__file__).resolve().parents[1]
    with zipfile.ZipFile(provenance / "source.zip","w",zipfile.ZIP_DEFLATED) as archive:
        for folder in ("src","research"):
            for path in (root / folder).rglob("*"):
                if path.is_file() and "runs" not in path.parts and "__pycache__" not in path.parts:
                    archive.write(path,path.relative_to(root))
        archive.write(root / "CMakeLists.txt","CMakeLists.txt")
    summary={"status":"running","arguments":{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
        "model":"f0=(1-epsilon)N(0,I)+epsilon N(0,diag(2,.5,.5)); M=N(0,I), g0=epsilon(B-M)",
        "resampling":False,"synchronization":False,"advection":False,"projection":False,
        "maxwellian":"fixed exact mass=1, mean=0, scalar temperature=1; collisions conserve these continuum invariants",
        "executable_sha256":hashlib.sha256(executable.read_bytes()).hexdigest(),"cases":[]}
    summary_path=args.output / "summary.json"
    try:
        # Exact same-seed replay and a zero-collision control validate the harness.
        control=dict(replicas=2,steps=2,full_count=args.full_count,sign_count=args.sign_count,
                     epsilon=args.epsilons[0],dt=args.dt,strength=0.,seed=41000)
        one=run(executable,args.output / "zero_control",**control)
        two=run(executable,args.output / "zero_replay",**control)
        if any(not np.array_equal(one['numerical'][key],two['numerical'][key]) for key in one['numerical']):
            raise RuntimeError("Same-seed replay failed")
        if any(not np.array_equal(one[key],np.repeat(one[key][:,:1,:],3,axis=1)) for key in ('full','signed')):
            raise RuntimeError("Zero-collision control changed its state")
        summary['zero_collision_and_replay_checks']='passed'
        for index,epsilon in enumerate(args.epsilons):
            directory=args.output / f"epsilon_{epsilon:g}"
            directory.mkdir()
            common=dict(steps=args.steps,full_count=args.full_count,sign_count=args.sign_count,
                        epsilon=epsilon,dt=args.dt,strength=args.strength)
            base=100000000+index*100000000
            print(f"epsilon={epsilon} pilot",flush=True)
            pilot=run(executable,directory / "pilot",replicas=args.pilot,seed=base,**common)
            print(f"epsilon={epsilon} evaluation",flush=True)
            evaluation=run(executable,directory / "evaluation",replicas=args.replicas,seed=base+10000000,**common)
            print(f"epsilon={epsilon} reference",flush=True)
            reference=run(executable,directory / "reference",replicas=args.reference_replicas,seed=base+20000000,
                          **{**common,'full_count':args.reference_count},mode="pic")
            print(f"epsilon={epsilon} refined reference",flush=True)
            refined=run(executable,directory / "reference_dt_half",replicas=args.reference_replicas,seed=base+30000000,
                        **{**common,'full_count':args.reference_count,'steps':2*args.steps,'dt':args.dt/2},mode="pic")
            weight,independent,vf,vd,cov=fit(pilot)
            f,d=evaluation['full'],evaluation['signed']
            start=time.perf_counter()
            mixed=weight*f+(1-weight)*d
            blend_time=time.perf_counter()-start
            independent_mix=independent*f+(1-independent)*d
            proxy=evaluation['proxy'][:,:,None]*f+(1-evaluation['proxy'][:,:,None])*d
            truth=reference['full'].mean(axis=0)
            truth_se=reference['full'].std(axis=0,ddof=1)/np.sqrt(args.reference_replicas)
            initial=np.array([1.5*epsilon,3+9*epsilon,epsilon*(np.exp(-1)-np.exp(-.25)),1,0,0,0,3.])
            truth[0]=initial; truth_se[0]=0.
            refined_mean=refined['full'][:,::2,:].mean(axis=0)
            refined_se=refined['full'][:,::2,:].std(axis=0,ddof=1)/np.sqrt(args.reference_replicas)
            methods={"full":f,"signed":d,"covariance_mix":mixed,"uncorrelated_mix":independent_mix,"count_proxy_mix":proxy}
            ci=bootstrap_gain(f[:,:,:3],d[:,:,:3],mixed[:,:,:3],reference['full'][:,:,:3],initial[:3],np.random.default_rng(base+40000000))
            records=[]
            for step in range(args.steps+1):
                for j,name in enumerate(NAMES):
                    statistics={method:{"mean":float(values[:,step,j].mean()),
                        "variance":float(values[:,step,j].var(ddof=1)),
                        "reference_bias":float(values[:,step,j].mean()-truth[step,j]),
                        "reference_mse":float(np.mean((values[:,step,j]-truth[step,j])**2))}
                        for method,values in methods.items()}
                    vfull=statistics['full']['variance']; vsigned=statistics['signed']['variance']
                    row={"step":step,"time":step*args.dt,"observable":name,"reference_mean":float(truth[step,j]),
                        "reference_mean_se":float(truth_se[step,j]),"reference_dt_half_difference":float(refined_mean[step,j]-truth[step,j]),
                        "reference_dt_difference_se":float(np.hypot(truth_se[step,j],refined_se[step,j])),
                        "pilot_full_weight":float(weight[step,j]),"pilot_independent_weight":float(independent[step,j]),
                        "pilot_covariance":float(cov[step,j]),
                        "evaluation_covariance":float(np.cov(f[:,step,j],d[:,step,j],ddof=1)[0,1]),
                        "evaluation_correlation":float(np.cov(f[:,step,j],d[:,step,j],ddof=1)[0,1]/np.sqrt(vfull*vsigned)) if vfull*vsigned>1e-25 else None,
                        "methods":statistics}
                    if j<3:
                        row.update(variance_gain=min(vfull,vsigned)/statistics['covariance_mix']['variance'],
                            reference_mse_gain=min(statistics['full']['reference_mse'],statistics['signed']['reference_mse'])/statistics['covariance_mix']['reference_mse'],
                            reference_mse_gain_ci95=[float(v) for v in ci[:,step,j]])
                    records.append(row)
            final_counts={name:{"mean":float(value[:,-1].mean()),"min":float(value[:,-1].min()),"max":float(value[:,-1].max())}
                          for name,value in evaluation['counts'].items()}
            case={"epsilon":epsilon,"records":records,"final_counts":final_counts,
                  "physical_cost_seconds_per_replica":evaluation['cost'],"blend_seconds_all_replicas":blend_time,
                  "pilot_wall_seconds":pilot['wall_seconds'],"evaluation_wall_seconds":evaluation['wall_seconds'],
                  "reference_wall_seconds":reference['wall_seconds'],"refined_reference_wall_seconds":refined['wall_seconds']}
            summary['cases'].append(case)
            summary_path.write_text(json.dumps(summary,indent=2)+"\n")
        summary['status']='complete'
    except Exception as error:
        summary.update(status='failed',failure=str(error))
        raise
    finally:
        summary_path.write_text(json.dumps(summary,indent=2)+"\n")
    report(summary,args.output)


def report(summary,output):
    lines=["# Homogeneous Coulomb mixing before any resampling", "",summary['model'],"",
        "Only the existing C++ homogeneous collisions and legacy signed source sampler evolve the state. No advection, projection, synchronization, adaptation or resampling is called. The Maxwellian has the exact continuum invariant moments and is fixed.","",
        "Independent pilot ensembles select an observable- and time-specific covariance-aware mixing weight. Fresh evaluation ensembles measure performance. Mixtures are diagnostics of the same evolved trajectories; mixing does not feed back into collisions.","",
        "References are independent high-count PIC ensembles. A separate half-timestep reference checks sensitivity; it is not an exact solution. Bootstrap intervals resample evaluation trajectories and reference trajectories, conditional on the fitted pilot weights.","",
        "Gain > 1 favors mixing. Comparisons reuse both representations already evolved by HDP, at the same physical cost.","",
        "| epsilon | Time | Observable | Variance gain | Reference MSE gain (95% interval) | Full weight | Full-signed correlation |",
        "|---:|---:|---|---:|---:|---:|---:|"]
    for case in summary['cases']:
        rows=case['records']; end=max(r['step'] for r in rows)
        for row in rows:
            if row['observable'] not in NAMES[:3] or row['step'] not in (0,end//2,end): continue
            low,high=row['reference_mse_gain_ci95']
            lines.append(f"| {case['epsilon']:g} | {row['time']:g} | {row['observable']} | {row['variance_gain']:.3f} | {row['reference_mse_gain']:.3f} [{low:.3f}, {high:.3f}] | {row['pilot_full_weight']:.3f} | {row['evaluation_correlation']:.3f} |")
    lines += ["", "## Final bias and reference sensitivity", "",
        "| epsilon | Observable | Full / signed / mixed reference bias | Reference mean SE | Half-dt reference change ± combined SE |",
        "|---:|---|---:|---:|---:|"]
    for case in summary['cases']:
        end=max(r['step'] for r in case['records'])
        for row in case['records']:
            if row['step']!=end or row['observable'] not in NAMES[:3]:continue
            biases='/'.join(f"{row['methods'][m]['reference_bias']:.4g}" for m in ('full','signed','covariance_mix'))
            lines.append(f"| {case['epsilon']:g} | {row['observable']} | {biases} | {row['reference_mean_se']:.3g} | {row['reference_dt_half_difference']:.3g} ± {row['reference_dt_difference_se']:.3g} |")
    lines += ["", "## Conservation, particle growth and timing", ""]
    for case in summary['cases']:
        lines.append(f"epsilon={case['epsilon']:g}: final count ranges `{case['final_counts']}`; per-replica initialization/collision/diagnostic costs `{case['physical_cost_seconds_per_replica']}`; blend over all evaluation replicas {case['blend_seconds_all_replicas']:.6g} seconds.")
        for name in NAMES[3:]:
            rows=[r for r in case['records'] if r['observable']==name]
            base=rows[0]
            shifts={m:max(abs(r['methods'][m]['mean']-base['methods'][m]['mean']) for r in rows) for m in ('full','signed','covariance_mix')}
            lines.append(f"- Maximum ensemble-mean change in {name}: `{shifts}`.")
    lines += ["", "The full reference MSE, variance, covariance, uncorrelated mixing and the existing count-proxy mixing results are in summary.json. All raw trajectories, seeds, executable and source provenance are retained. Pilot simulation costs are reported separately; amortizing those costs is required before claiming operational efficiency. This is an estimator test before resampling, not a validation of the full adaptive solver."]
    lines += ["", "## Final comparisons with simpler weights", "",
        "These are reference-MSE gains versus the better component. Only the covariance-aware weights use covariance in the pilot; the count proxy is computed from the evaluation population counts, as in the solver.", "",
        "| epsilon | Observable | Covariance-aware | Uncorrelated pilot | Count proxy |", "|---:|---|---:|---:|---:|"]
    for case in summary['cases']:
        end=max(r['step'] for r in case['records'])
        for row in case['records']:
            if row['step']!=end or row['observable'] not in NAMES[:3]: continue
            best=min(row['methods'][m]['reference_mse'] for m in ('full','signed'))
            gains=[best/row['methods'][m]['reference_mse'] for m in ('covariance_mix','uncorrelated_mix','count_proxy_mix')]
            lines.append(f"| {case['epsilon']:g} | {row['observable']} | {gains[0]:.3f} | {gains[1]:.3f} | {gains[2]:.3f} |")
    (output / "REPORT.md").write_text("\n".join(lines)+"\n",encoding='utf-8')


if __name__=="__main__":main()
