"""Measure homogeneous HDP/PIC cost ratios at a common relative-error threshold."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile
import numpy as np
from matched_accuracy import metrics, read_run, reference_estimates, run, scales


def combine(batches):
    return dict(settings=batches[0]['settings'],
                estimates={key: np.concatenate([b['estimates'][key] for b in batches])
                           for key in batches[0]['estimates']},
                costs=np.concatenate([b['costs'] for b in batches]))


def ratio_interval(pic_batches, hdp_batches, seed, bootstraps=2000):
    """Resample sequential timing batches, then trajectories within batches."""
    rng = np.random.default_rng(seed)
    def draw(batches):
        selected = rng.integers(len(batches), size=len(batches))
        return np.mean([rng.choice(batches[i], size=len(batches[i])).mean() for i in selected])
    ratios = [draw(hdp_batches) / draw(pic_batches) for _ in range(bootstraps)]
    return np.quantile(ratios, [.025, .975]).tolist()


def plot(summary, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = [c for c in summary['cases'] if 'ratio' in c]
    if not rows:
        return
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
    for row in rows:
        x, y = row['epsilon'], row['ratio']
        lo, hi = row['ratio_ci95']
        axes[0].errorbar(x, y, yerr=[[max(0, y-lo)], [max(0, hi-y)]],
                        fmt='o' if row['verified'] else 'x', color='#1565a5' if row['verified'] else '#c23b32',
                        capsize=4, markersize=7)
        for method, color, marker in (('pic', '#444444', 's'), ('hdp', '#1565a5', 'o')):
            m = row[method]['metrics']
            v = m['joint_relative_rmse']; low, high = m['joint_relative_rmse_ci95']
            axes[1].errorbar(x, v, yerr=[[max(0, v-low)], [max(0, high-v)]],
                            fmt=marker, color=color, capsize=3)
    verified = [r for r in rows if r['verified']]
    if verified:
        x = np.array([r['epsilon'] for r in verified]); y = np.array([r['ratio'] for r in verified])
        axes[0].plot(x, y, color='#1565a5', alpha=.55, label='Validated common threshold')
        anchor = len(x)//2
        axes[0].plot(x, y[anchor]*x/x[anchor], '--', color='#ac7b20', label='Proportional to epsilon (guide)')
    if any(not r['verified'] for r in rows):
        axes[0].plot([], [], 'x', color='#c23b32', label='Accuracy target not verified')
    axes[0].axhline(1, color='black', linestyle=':', label='Equal cost')
    axes[0].set_yscale('log')
    axes[0].set_ylabel('R = HDP + mixing time / PIC time')
    axes[0].set_title('Homogeneous collisions: ordinary estimators\nR < 1 favors HDP + mixing')
    axes[0].legend(fontsize=9)
    axes[1].axhline(summary['arguments']['target'], color='black', linestyle=':', label='Common accuracy threshold')
    axes[1].plot([], [], 's', color='#444444', label='PIC')
    axes[1].plot([], [], 'o', color='#1565a5', label='HDP + mixing')
    axes[1].set_ylabel('Normalized trajectory RMS error')
    axes[1].set_xlabel('Perturbation size epsilon (log scale)')
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.set_xscale('log'); ax.grid(alpha=.2, which='both')
    axes[1].set_xticks([r['epsilon'] for r in rows], [f"{r['epsilon']:.3g}" for r in rows])
    fig.savefig(output/'R_epsilon.png', dpi=200)
    fig.savefig(output/'R_epsilon.svg')
    plt.close(fig)


def report(summary, output):
    lines = ['# Homogeneous efficiency as a function of perturbation size', '',
        'R(epsilon) = complete HDP-plus-mixing compute time / ordinary PIC compute time. Lower values favor HDP. Both methods use the same three-observable normalized trajectory RMS target. Particle allocations are selected separately using pilot upper bootstrap bounds at 90% of the final threshold, then tested on fresh seeds. This measures a common accuracy threshold, not exact equality of achieved errors or globally optimal counts.', '',
        'T=0.2, dt=0.01, collision strength 5. No resampling, projection, spatial evolution or invariant correction. Compute costs include initialization, source bounds, collisions, signed-source sampling and diagnostics/mixing; CSV output is excluded. All competitors run sequentially. Three timing batches alternate method order. Ratio intervals hierarchically bootstrap batches and trajectories; only three batches are available, so these intervals cannot capture every machine-load effect.', '',
        '| epsilon | PIC / HDP RMS | PIC count | HDP full / per sign | PIC / HDP time (s) | R [95% interval] | Accuracy verified |',
        '|---:|---|---:|---|---|---|---|']
    for c in summary['cases']:
        if 'ratio' not in c:
            lines.append(f"| {c['epsilon']:.5g} | unavailable | | | | | False |")
            continue
        p, h = c['pic'], c['hdp']; pm, hm = p['metrics'], h['metrics']
        lo, hi = c['ratio_ci95']
        lines.append(f"| {c['epsilon']:.5g} | {pm['joint_relative_rmse']:.3f} / {hm['joint_relative_rmse']:.3f} | {p['settings']['full_count']} | {h['settings']['full_count']} / {h['settings']['sign_count']} | {pm['mean_compute_seconds']:.5g} / {hm['mean_compute_seconds']:.5g} | {c['ratio']:.3f} [{lo:.3f}, {hi:.3f}] | {c['verified']} |")
    lines += ['', 'Reference uncertainty and doubled-count/half-timestep differences, all pilot allocations, raw validation batches and failed accuracy checks are retained in summary.json and run directories. Intervals are conditional on selected allocations and approximate references. Kernel approximation, finite radial support and empirical rejection envelopes remain limitations. The proportional-to-epsilon line in the figure is a visual guide, not a fitted law or proof of asymptotic scaling.']
    (output/'REPORT.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    plot(summary, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--epsilon-min', type=float, default=.01)
    parser.add_argument('--epsilon-max', type=float, default=.3)
    parser.add_argument('--points', type=int, default=7)
    parser.add_argument('--target', type=float, default=.4)
    parser.add_argument('--pilot-replicas', type=int, default=48)
    parser.add_argument('--batch-replicas', type=int, default=64)
    args = parser.parse_args()
    if not (0 < args.epsilon_min < args.epsilon_max < 1) or args.points < 2 or args.target <= 0:
        parser.error('Require 0 < epsilon-min < epsilon-max < 1, points >= 2, target > 0')
    if min(args.pilot_replicas, args.batch_replicas) < 8:
        parser.error('Use at least eight trajectories per batch')
    output = args.output
    arguments = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items() if k != 'resume'}
    if args.resume:
        previous = json.loads((output/'summary.json').read_text())
        if previous['arguments'] != arguments:
            raise ValueError('Resume arguments differ from recorded protocol')
    else:
        output.mkdir(parents=True)
        provenance = output/'provenance'; provenance.mkdir()
        shutil.copy2(args.executable, provenance/args.executable.name)
        for dll in args.executable.parent.glob('*.dll'):
            shutil.copy2(dll, provenance/dll.name)
        root = Path(__file__).resolve().parents[1]
        with zipfile.ZipFile(provenance/'source.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
            for folder in ('src', 'research', 'tests'):
                for path in (root/folder).rglob('*'):
                    if path.is_file() and 'runs' not in path.parts and '__pycache__' not in path.parts:
                        archive.write(path, path.relative_to(root))
            archive.write(root/'CMakeLists.txt', 'CMakeLists.txt')
    executable = (output/'provenance'/args.executable.name).resolve()
    summary = dict(status='running', arguments=arguments, cases=[],
                   executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
    def save():
        (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    save()
    try:
        for ci, epsilon in enumerate(np.geomspace(args.epsilon_min, args.epsilon_max, args.points)):
            epsilon = float(epsilon)
            folder = output/f'epsilon_{ci}_{epsilon:.8g}'; folder.mkdir(exist_ok=True)
            base = 150000000 + ci*350000000
            def config(mode, nf, ns, replicas, offset, dt=.01, steps=20):
                return dict(mode=mode, full_count=nf, sign_count=ns, replicas=replicas,
                            seed=base+offset, epsilon=epsilon, dt=dt, steps=steps, strength=5.)
            print(f'epsilon={epsilon:.6g}: references', flush=True)
            reference = reference_estimates(run(executable, folder/'reference', config('pic',32768,16,128,0)), epsilon)
            double = reference_estimates(run(executable, folder/'reference_double', config('pic',65536,16,32,2000000)),epsilon)
            half = reference_estimates(run(executable, folder/'reference_half_dt', config('pic',32768,16,32,4000000,.005,40)),epsilon)
            norm = scales(epsilon); truth = reference.mean(axis=0)
            case = dict(epsilon=epsilon, candidates=[],
                reference_uncertainty=float(np.sqrt(np.mean((reference.std(axis=0,ddof=1)/np.sqrt(len(reference))/norm)**2))),
                reference_count_difference=float(np.sqrt(np.mean(((double.mean(axis=0)-truth)/norm)**2))),
                reference_dt_difference=float(np.sqrt(np.mean(((half[:,::2].mean(axis=0)-truth)/norm)**2))))
            summary['cases'].append(case)
            def candidate(mode, nf, ns):
                label = f'pilot_{mode}_f{nf}_s{ns}'
                print(f'epsilon={epsilon:.6g}: {label}', flush=True)
                data = run(executable, folder/label, config(mode,nf,ns,args.pilot_replicas,10000000+len(case['candidates'])*2000000))
                method = 'full' if mode == 'pic' else 'mixed'
                m = metrics(data,reference,epsilon,base+70000000+len(case['candidates']))[method]
                row = dict(label=label,settings=data['settings'],metrics=m)
                case['candidates'].append(row); save()
            selection_target = .9*args.target
            for ns in (32,64,128,256):
                for ratio in (8,16):
                    candidate('hdp',ns*ratio,ns)
            # The inverse-square estimate only proposes counts; all points are measured.
            estimate = 2/(epsilon*epsilon*selection_target*selection_target)
            for nf in sorted({max(128, int(np.ceil(estimate*f/32))*32) for f in (.6,1.,1.6)}):
                candidate('pic',nf,16)
            selected = {}
            for mode in ('pic','hdp'):
                feasible = [r for r in case['candidates'] if r['settings']['mode']==mode and
                            r['metrics']['joint_relative_rmse_ci95'][1] <= selection_target]
                if not feasible and mode=='pic':
                    nf = max(r['settings']['full_count'] for r in case['candidates'] if r['settings']['mode']=='pic')*2
                    candidate('pic',nf,16)
                    feasible = [r for r in case['candidates'] if r['settings']['mode']==mode and r['metrics']['joint_relative_rmse_ci95'][1]<=selection_target]
                if feasible:
                    selected[mode] = min(feasible,key=lambda r:r['metrics']['mean_compute_seconds'])
            if len(selected)<2:
                case['failure']='No feasible pilot allocation in measured search'; save(); report(summary,output); continue
            batches = {'pic':[], 'hdp':[]}
            for batch in range(3):
                for mode in (('pic','hdp') if batch%2==0 else ('hdp','pic')):
                    chosen = selected[mode]['settings']
                    offset = 100000000+batch*10000000+(2000000 if mode=='hdp' else 0)
                    label=f'validation_{mode}_batch{batch}'
                    print(f'epsilon={epsilon:.6g}: {label}',flush=True)
                    data=run(executable,folder/label,config(mode,chosen['full_count'],chosen['sign_count'],args.batch_replicas,offset))
                    batches[mode].append(data)
            for mode in ('pic','hdp'):
                merged=combine(batches[mode]); method='full' if mode=='pic' else 'mixed'
                case[mode]=dict(settings=selected[mode]['settings'],validation_replicas=len(merged['costs']),
                    metrics=metrics(merged,reference,epsilon,base+190000000)[method],
                    batch_mean_seconds=[float(b['costs'].mean()) for b in batches[mode]])
            case['verified']=all(case[mode]['metrics']['joint_relative_rmse_ci95'][1]<=args.target for mode in ('pic','hdp'))
            case['ratio']=case['hdp']['metrics']['mean_compute_seconds']/case['pic']['metrics']['mean_compute_seconds']
            case['ratio_ci95']=ratio_interval([b['costs'] for b in batches['pic']],[b['costs'] for b in batches['hdp']],base+195000000)
            save(); report(summary,output)
        summary['status']='complete'
    except Exception as error:
        summary.update(status='failed',failure=str(error)); raise
    finally:
        save()
    report(summary,output)
    print(output/'R_epsilon.png',flush=True)


if __name__=='__main__':
    main()
