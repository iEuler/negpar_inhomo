"""Fresh conservative confirmations and paired mixing gains for the repaired sweep."""
import argparse
import json
from pathlib import Path
import numpy as np
from matched_accuracy import metrics, read_run, reference_estimates, run, scales


def mixing_gain(data, reference, epsilon, seed, bootstraps=1000):
    """Compare the mixture to components of the same evolving populations."""
    names = ('full', 'signed', 'mixed')
    values = np.stack([data['estimates'][name] for name in names]) / scales(epsilon)
    reference = reference / scales(epsilon)
    truth = reference.mean(axis=0)
    mse = ((values - truth) ** 2).mean(axis=(1, 2, 3))
    rng = np.random.default_rng(seed)
    gains = []
    for _ in range(bootstraps):
        indices = rng.integers(values.shape[1], size=values.shape[1])
        ref = reference[rng.integers(len(reference), size=len(reference))].mean(axis=0)
        losses = ((values[:, indices] - ref) ** 2).mean(axis=(1, 2, 3))
        gains.append([losses[1] / losses[2], min(losses[:2]) / losses[2]])
    bounds = np.quantile(gains, [.025, .975], axis=0)
    return dict(component_relative_rmse=dict(zip(names, map(float, np.sqrt(mse)))),
                signed_over_mixed_mse=float(mse[1] / mse[2]),
                signed_over_mixed_ci95=bounds[:, 0].tolist(),
                better_component_over_mixed_mse=float(min(mse[:2]) / mse[2]),
                better_component_over_mixed_ci95=bounds[:, 1].tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('experiment', type=Path)
    args = parser.parse_args()
    out = args.experiment
    sweep = json.loads((out / 'summary.json').read_text())
    if sweep['status'] != 'complete':
        raise ValueError('Complete the sweep before supplemental analysis')
    confirmation = out / 'confirmation'
    confirmation.mkdir(exist_ok=True)
    executable = (out / 'provenance/negpar_homogeneous.exe').resolve()
    # Allocations fixed after the exploratory sweep, before these fresh runs.
    # Keep its failed validations intact; these are supplemental comparisons.
    allocations = [(.05, 'ordinary', 8192, 2048, 128),
                   (.05, 'initial-moment CV', 1024, 256, 32),
                   (.3, 'ordinary closer error', 768, 512, 32)]
    result = dict(confirmations=[], mixing=[])
    for index, (epsilon, label, nf_pic, nf_hdp, np_hdp) in enumerate(allocations):
        folder = out / f'epsilon_{epsilon:g}'
        reference = reference_estimates(read_run(folder / 'reference'), epsilon)
        measured = []
        settings = []
        for side, (mode, nf, np_) in enumerate((('pic', nf_pic, 16), ('hdp', nf_hdp, np_hdp))):
            config = dict(mode=mode, full_count=nf, sign_count=np_, replicas=384,
                          seed=2400000000 + index * 10000000 + side * 2000000,
                          epsilon=epsilon, dt=.01, steps=20, strength=5.)
            directory = confirmation / f'epsilon{epsilon:g}_{index}_{mode}'
            print(directory.name, flush=True)
            data = run(executable, directory, config)
            method = ('pic_cv' if mode == 'pic' else 'mixed_cv') if label == 'initial-moment CV' else ('full' if mode == 'pic' else 'mixed')
            measured.append(metrics(data, reference, epsilon, config['seed'] + 77)[method])
            settings.append(config)
        p, h = measured
        result['confirmations'].append(dict(epsilon=epsilon, comparison=label, target=.4,
            pic=p, hdp=h, pic_settings=settings[0], hdp_settings=settings[1],
            target_verified=all(x['joint_relative_rmse_ci95'][1] <= .4 for x in measured),
            time_ratio=p['mean_compute_seconds'] / h['mean_compute_seconds']))
        (confirmation / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    # Include all sweep candidates and fresh validations. Individual intervals
    # are descriptive, not simultaneous tests or a selection-corrected claim.
    for ci, case in enumerate(sweep['cases']):
        folder = out / f"epsilon_{case['epsilon']:g}"
        reference = reference_estimates(read_run(folder / 'reference'), case['epsilon'])
        directories = sorted(folder.glob('hdp_*')) + sorted(folder.glob('validate_hdp_*'))
        directories += sorted(confirmation.glob(f"epsilon{case['epsilon']:g}_*_hdp"))
        for di, directory in enumerate(directories):
            data = read_run(directory)
            row = dict(epsilon=case['epsilon'], directory=str(directory), settings=data['settings'],
                       gain=mixing_gain(data, reference, case['epsilon'], 170000 + ci * 1000 + di))
            raw = np.genfromtxt(directory / 'trajectories.csv', delimiter=',', names=True)
            row['signed_invariant_drift'] = {}
            for name, j in (('mass', 3), ('px', 4), ('py', 5), ('pz', 6), ('v2', 7)):
                values = raw[f'signed_{j}'].reshape(data['settings']['replicas'], -1)
                change = values[:, -1] - values[:, 0]
                row['signed_invariant_drift'][name] = dict(mean=float(change.mean()),
                    se=float(change.std(ddof=1) / np.sqrt(len(change))))
            result['mixing'].append(row)
    (confirmation / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = ['# Fresh confirmations after the repaired-source sweep', '',
        'Allocations were fixed after the exploratory sweep and evaluated with 384 new seeds per method. Original failed validations remain in the main report. Costs include all homogeneous solver work and exclude CSV output. These are supplemental measured allocations, not a new optimization.', '',
        '| epsilon | Comparison | PIC RMS [95% interval] | HDP RMS [95% interval] | PIC / HDP seconds | Time ratio | Target verified |',
        '|---:|---|---|---|---|---:|---|']
    for row in result['confirmations']:
        p, h = row['pic'], row['hdp']
        fmt = lambda x: f"{x['joint_relative_rmse']:.3f} [{x['joint_relative_rmse_ci95'][0]:.3f}, {x['joint_relative_rmse_ci95'][1]:.3f}]"
        lines.append(f"| {row['epsilon']:g} | {row['comparison']} | {fmt(p)} | {fmt(h)} | {p['mean_compute_seconds']:.5f} / {h['mean_compute_seconds']:.5f} | {row['time_ratio']:.3f} | {row['target_verified']} |")
    lines += ['', 'Time ratio > 1 favors HDP. The common threshold is 0.4; achieved errors differ. Reference uncertainty, count and timestep checks are in the main report. Small runtime differences need repeated timing measurements before a strong efficiency claim.', '',
        '# Incremental mixing on the same populations', '',
        'Gain is normalized trajectory MSE of a component divided by mixture MSE. Paired bootstrap resamples whole HDP trajectories and independent PIC reference trajectories (1,000 replicates). The better component is chosen on aggregate normalized MSE. Costs are shared because the components already exist; this does not measure a standalone signed-only solver. Intervals are individual, without multiple-comparison or selection correction.', '',
        '| epsilon | Run | Signed / mixture MSE [95% interval] | Better component / mixture MSE [95% interval] |',
        '|---:|---|---|---|']
    for row in result['mixing']:
        g = row['gain']
        fmt = lambda key: f"{g[key + '_mse']:.3f} [{g[key + '_ci95'][0]:.3f}, {g[key + '_ci95'][1]:.3f}]"
        lines.append(f"| {row['epsilon']:g} | {Path(row['directory']).name} | {fmt('signed_over_mixed')} | {fmt('better_component_over_mixed')} |")
    lines += ['', 'All signed invariant drift means and standard errors are saved in summary.json. Passing rejection checks does not certify the empirical envelope globally, finite-support error, approximate kernel, or exact expected conservation.']
    (confirmation / 'REPORT.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    (out / 'provenance/split_sweep_analysis.py').write_bytes(Path(__file__).read_bytes())
    print(confirmation / 'REPORT.md')


if __name__ == '__main__':
    main()
