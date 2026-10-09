"""Summarize the source audit and clearly supersede the affected old runs."""
import json
from pathlib import Path
import numpy as np
from matched_accuracy import read_run, reference_estimates, initial
from homogeneous_experiment import fit, bootstrap_gain


def main():
    runs=Path(__file__).parent/'runs'
    notice=('**Superseded: source-cache bug.** These HDP trajectories used zero cached '
            'signed counts in the source rejection denominator. Their conservation '
            'and efficiency conclusions are withdrawn pending corrected measurements. '
            'Pure PIC trajectories and frozen-distribution tests are unaffected. '
            'Raw files are preserved. See `../source_conservation_corrected_v1/REPORT.md`.')
    for label in ('homogeneous_v1','homogeneous_dt_half_v1','matched_accuracy_v1'):
        directory=runs/label
        if not directory.exists():continue
        (directory/'VALIDITY.md').write_text(notice+'\n')
        for report in directory.glob('*REPORT.md'):
            text=report.read_text(encoding='utf-8')
            if '**Superseded: source-cache bug.**' not in text:report.write_text(notice+'\n\n'+text,encoding='utf-8')
        if label=='matched_accuracy_v1':
            report=directory/'confirmation/REPORT.md'
            if report.exists():
                text=report.read_text()
                if '**Superseded: source-cache bug.**' not in text:
                    report.write_text(notice.replace('../source_conservation','../../source_conservation')+'\n\n'+text)
    out=runs/'source_conservation_corrected_v1'
    summary=json.loads((out/'summary.json').read_text())
    stages=json.loads((runs/'source_conservation_v1/summary.json').read_text())['stages']
    lines=['# Signed-source conservation audit', '',
           'Root cause: `samplefromhNeg` computed its rejection denominator from cached `positiveMoments.m0` and `negativeMoments.m0`. `pushBack` and merging do not update these caches. The homogeneous harness never computed them, so both were zero. Nonzero targets divided by zero and were effectively accepted with probability one; zero targets produced NaN and were rejected. The proposal count still used live sizes. This changes the sampled distribution, not just its variance.', '',
           'The production fix uses live particle-list sizes for both rejection densities. No source formula, collision kernel, moment projection or resampling was changed. A seeded regression requires identical samples with fresh and stale caches.', '',
           '## Paired one-step isolation (2,048 replicas per timestep)', '',
           '| dt | Cache | Source mass change (mean ± SE) | Source plus transport total-v² change (mean ± SE) |',
           '|---:|---|---:|---:|']
    for row in stages:
        fmt=lambda key:f"{row[key]['mean']:.6g} ± {row[key]['se']:.3g}"
        lines.append(f"| {row['dt']} | {'fresh' if row['fresh_cache'] else 'stale'} | {fmt('source_mass')} | {fmt('combined_v2')} |")
    lines += ['', 'Transport mass change is exactly zero in all replicas: that stage changes velocities, not signed counts. Source and transport energy changes must be assessed together; neither conserves signed energy separately.', '',
              '## Corrected trajectories, final time 0.2', '',
              '| Configuration | Replicas | Signed mass drift (mean ± SE) | Signed total-v² drift (mean ± SE) |',
              '|---|---:|---:|---:|']
    for row in summary:
        fmt=lambda key:f"{row['drift'][key]['mean']:.6g} ± {row['drift'][key]['se']:.3g}"
        lines.append(f"| {row['label']} | {row['settings']['replicas']} | {fmt('signed_mass')} | {fmt('signed_total_v2')} |")
    lines += ['', 'Mass and total-v² ensemble drifts are consistent with zero at this precision. One momentum component in the epsilon=0.3 half-dt run is -0.00171 ± 0.000585 SE (about 2.9 SE from zero); multiple component checks and finite samples do not certify all invariant expectations. All component measurements are retained in summary.json.', '',
              '## Selected efficiency points after the fix', '',
              'PIC measurements are reused from the unaffected frozen pure-PIC runs. HDP measurements use a new frozen corrected executable and fresh seeds. These are selected-point rechecks, not a new count optimization.', '',
              '| Estimator at epsilon=0.01 | PIC RMS | Corrected HDP RMS (95% interval) | PIC / HDP compute time |',
              '|---|---:|---:|---:|']
    for row in summary:
        if 'pic_comparison' not in row:continue
        pair=row['pic_comparison'];h=row['metrics'][pair['hdp_method']];p=pair['pic']
        ci=h['joint_relative_rmse_ci95']
        lines.append(f"| {pair['hdp_method']} | {p['joint_relative_rmse']:.3f} | {h['joint_relative_rmse']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] | {pair['time_ratio']:.2f} |")
    lines += ['', '## Mixing recheck using corrected trajectories', '',
              'The first 256 trajectories fit weights; the remaining 768 evaluate them. Pure-PIC references use known initial expectations as a temporal control variate. Gains compare final-time MSE to the better component, reusing the same populations. Intervals resample evaluation and reference trajectories, conditional on pilot weights.', '',
              '| epsilon | Observable | Covariance-mixture gain (95% interval) | Count-proxy gain (95% interval) |',
              '|---:|---|---:|---:|']
    rechecks=[]
    for epsilon,label in ((.05,'eps0.05'),(.3,'eps0.3')):
        data=read_run(out/label)
        f,d=data['estimates']['full'],data['estimates']['signed']
        weights=fit(dict(full=f[:256],signed=d[:256]))[0]
        f,d=f[256:],d[256:]
        mixed=weights*f+(1-weights)*d
        proxy=data['estimates']['mixed'][256:]
        ref=reference_estimates(read_run(runs/f'matched_accuracy_v1/epsilon_{epsilon:g}/reference'),epsilon)
        truth=ref.mean(axis=0)
        for name,values in (('covariance',mixed),('proxy',proxy)):
            ci=bootstrap_gain(f,d,values,ref,initial(epsilon),np.random.default_rng(2026000+int(epsilon*10000)))
            mse=lambda x:np.mean((x[:,-1]-truth[-1])**2,axis=0)
            gain=np.minimum(mse(f),mse(d))/mse(values)
            rechecks.append(dict(epsilon=epsilon,method=name,gain=gain.tolist(),ci95=ci[:,-1,:].tolist()))
        covariance,proxy_row=rechecks[-2:]
        for j,name in enumerate(('anisotropy','fourth moment','cosine difference')):
            fmt=lambda r:f"{r['gain'][j]:.3f} [{r['ci95'][0][j]:.3f}, {r['ci95'][1][j]:.3f}]"
            lines.append(f"| {epsilon} | {name} | {fmt(covariance)} | {fmt(proxy_row)} |")
    (out/'mixing_recheck.json').write_text(json.dumps(rechecks,indent=2)+'\n')
    lines += ['', '## Remaining source-sampler issues', '',
              'Source-centered spherical quadrature resolves the Coulomb singularity much better than Gaussian product quadrature. At radial order 256, the legacy kernel mass integrals are approximately -3.7e-7 and -1.6e-6 for source speeds 1 and 2. These small approximation errors do not explain the previous percent-level drift.', '',
              'A separate algebraic discrepancy remains in the legacy split. Write q=h-M and a=alphaNeg*M. The low branch includes q only where q<a, while the high branch samples max(q-a,0). Their sum is q-a*1(q>=a), rather than q. This omits a positive slice. Source-centered quadrature estimates omitted mass around 0.000684 and 0.000634 per unit source for speeds 1 and 2 at dt=0.01, before possible rejection-envelope clipping. Its effect on a positive-minus-negative source depends on both source populations; the corrected conservation tests quantify the net error without claiming global consistency.', '',
              'Empirical rejection bounds, finite radial support and the density used in positive proposal counts also require a separate envelope/normalization audit. They have not been altered here. Count dependence, long-time behavior and convergence remain unverified. The cache fix establishes neither exact trajectory conservation nor a globally conservative source sampler.', '',
              'Frozen binaries, sources and raw data are preserved in this directory and `../source_conservation_v1`. Tests: all 77 numerical cases passed, with 4,198 assertions, including the new cache-independence regression.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(out/'REPORT.md')


if __name__=='__main__':main()
