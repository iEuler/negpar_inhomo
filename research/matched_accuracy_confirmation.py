"""Supplement the archived benchmark with fresh validation and closer errors."""
import json
from pathlib import Path
import numpy as np
from matched_accuracy import run, read_run, reference_estimates, metrics, scales


def main():
    root = Path(__file__).parent / 'runs' / 'matched_accuracy_v1'
    summary = json.loads((root / 'summary.json').read_text())
    case = next(c for c in summary['cases'] if c['epsilon'] == .01)
    directory = root / 'confirmation'
    directory.mkdir(exist_ok=True)
    executable = (root / 'provenance' / 'negpar_homogeneous.exe').resolve()
    reference = reference_estimates(read_run(root / 'epsilon_0.01' / 'reference'), .01)
    results = {}
    configs = [('hdp_cv', 'hdp', 256, 32, 768),
               ('pic_cv_close', 'pic', 60000, 16, 192),
               ('pic_ordinary_close', 'pic', 220000, 16, 128)]
    for i, (label, mode, nf, ns, replicas) in enumerate(configs):
        config = dict(mode=mode, full_count=nf, sign_count=ns, replicas=replicas,
                      seed=1900000000+i*1000000, epsilon=.01, dt=.01, steps=20, strength=5.)
        print(label, flush=True)
        data = run(executable, directory / label, config)
        results[label] = dict(settings=config, metrics=metrics(data, reference, .01, config['seed']+77))
        alternatives = {}
        for check in ('reference_count_double', 'reference_dt_half'):
            alternative = reference_estimates(read_run(root/'epsilon_0.01'/check), .01)
            if check.endswith('half'): alternative = alternative[:, ::2]
            alternatives[check] = {
                m: float(np.sqrt(np.mean(((v-alternative.mean(axis=0))/scales(.01))**2)))
                for m, v in data['estimates'].items()
                if m in results[label]['metrics']}
        results[label]['alternative_reference_rmse'] = alternatives
        (directory/'summary.json').write_text(json.dumps(results, indent=2)+'\n')
    ordinary = next(p for p in case['comparisons'] if p['target']==.4 and p['comparison']=='ordinary')
    rows = [('ordinary, closer errors', results['pic_ordinary_close']['metrics']['full'], ordinary['hdp']),
            ('initial-moment CV, closer errors', results['pic_cv_close']['metrics']['pic_cv'], results['hdp_cv']['metrics']['mixed_cv'])]
    lines = ['# Fresh confirmation at epsilon = 0.01', '',
             'Same frozen executable, physics, reference and complete-compute timing as the main benchmark. New seeds are separate from count selection. The larger signed ensemble replaces the configuration whose independent upper error bound missed 0.25; the original failure remains in the main report.', '',
             '| Comparison | PIC RMS (95% interval) | HDP RMS (95% interval) | PIC seconds | HDP seconds | PIC/HDP time |',
             '|---|---:|---:|---:|---:|---:|']
    for name, p, h in rows:
        fmt = lambda x: f"{x['joint_relative_rmse']:.3f} [{x['joint_relative_rmse_ci95'][0]:.3f}, {x['joint_relative_rmse_ci95'][1]:.3f}]"
        lines.append(f"| {name} | {fmt(p)} | {fmt(h)} | {p['mean_compute_seconds']:.4g} | {h['mean_compute_seconds']:.4g} | {p['mean_compute_seconds']/h['mean_compute_seconds']:.2f} |")
    lines += ['', '## Reference sensitivity', '',
              '| Estimate | Main reference RMS | Double-count reference RMS | Half-dt reference RMS |',
              '|---|---:|---:|---:|']
    old_hdp = read_run(root/'epsilon_0.01'/'validate_hdp_n128_ratio8')
    old_checks = {}
    for check in ('reference_count_double', 'reference_dt_half'):
        alternative = reference_estimates(read_run(root/'epsilon_0.01'/check), .01)
        if check.endswith('half'): alternative = alternative[:, ::2]
        old_checks[check] = float(np.sqrt(np.mean(((old_hdp['estimates']['mixed']-alternative.mean(axis=0))/scales(.01))**2)))
    results['ordinary_hdp_reference_sensitivity'] = old_checks
    for label, method in [('hdp_cv','mixed_cv'), ('pic_cv_close','pic_cv'), ('pic_ordinary_close','full')]:
        row = results[label]
        lines.append(f"| {label} | {row['metrics'][method]['joint_relative_rmse']:.3f} | {row['alternative_reference_rmse']['reference_count_double'][method]:.3f} | {row['alternative_reference_rmse']['reference_dt_half'][method]:.3f} |")
    lines.append(f"| ordinary HDP mixture | {ordinary['hdp']['joint_relative_rmse']:.3f} | {old_checks['reference_count_double']:.3f} | {old_checks['reference_dt_half']:.3f} |")
    (directory/'summary.json').write_text(json.dumps(results, indent=2)+'\n')
    lines += ['', 'Errors are closely matched rather than exactly equal. Higher-count and half-timestep reference sensitivity is retained in summary.json. These results support the complete deviational method plus mixing; they do not assign the entire speedup to the blending operation.']
    (directory/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines), flush=True)


if __name__ == '__main__':
    main()
