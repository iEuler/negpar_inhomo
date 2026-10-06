"""Controlled mixing test with an exact frozen Gaussian-mixture distribution.

Requires NumPy; plot with frozen_mixing_plot.py. This tests scalar estimator
mixing, not the evolving C++ solver or its count-based mixing proxy.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np

OBSERVABLES = ("velocity", "kinetic_second_moment", "cosine_mode")


def moments(epsilon, shift, frequency):
    """Exact observable means and per-draw variances for f, G_shift, G_0."""
    def gaussian(mu):
        c = np.exp(-frequency**2 / 2) * np.cos(frequency * mu)
        return np.array([mu, 1 + mu**2, c]), np.array([
            1., 2 + 4 * mu**2,
            .5 * (1 + np.exp(-2 * frequency**2) * np.cos(2 * frequency * mu)) - c*c])
    m0, v0 = gaussian(0)
    m1, v1 = gaussian(shift)
    truth = (1-epsilon)*m0 + epsilon*m1
    vf = (1-epsilon)*(v0 + m0*m0) + epsilon*(v1 + m1*m1) - truth*truth
    return truth, m0, vf, epsilon**2 * (v0 + v1)


def evaluate(values, frequency):
    return np.column_stack((values.mean(axis=1), (values*values).mean(axis=1),
                            np.cos(frequency*values).mean(axis=1)))


def estimates(repetitions, full_count, sign_count, epsilon, shift, frequency, seeds, batch):
    full_rng, positive_rng, negative_rng = [np.random.default_rng(s) for s in seeds]
    full, signed = [], []
    full_time = signed_time = 0.
    baseline = moments(epsilon, shift, frequency)[1]
    for begin in range(0, repetitions, batch):
        n = min(batch, repetitions-begin)
        start = time.perf_counter()
        values = full_rng.standard_normal((n, full_count))
        values += shift * (full_rng.random(values.shape) < epsilon)
        full.append(evaluate(values, frequency))
        full_time += time.perf_counter()-start
        start = time.perf_counter()
        positive = positive_rng.standard_normal((n, sign_count)) + shift
        negative = negative_rng.standard_normal((n, sign_count))
        signed.append(baseline + epsilon*(evaluate(positive, frequency)-evaluate(negative, frequency)))
        signed_time += time.perf_counter()-start
    return np.concatenate(full), np.concatenate(signed), full_time, signed_time


def stats(values, truth):
    return {"mean": float(values.mean()), "bias": float(values.mean()-truth),
            "mean_standard_error": float(values.std(ddof=1)/np.sqrt(len(values))),
            "variance": float(values.var(ddof=1)), "mse": float(np.mean((values-truth)**2))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=20000)
    parser.add_argument("--pilot", type=int, default=4000)
    parser.add_argument("--full-count", type=int, default=512)
    parser.add_argument("--sign-count", type=int, default=256)
    parser.add_argument("--seed", type=int, default=84017)
    parser.add_argument("--batch", type=int, default=256)
    args = parser.parse_args()
    if min(args.repetitions, args.pilot) < 100 or min(args.full_count, args.sign_count, args.batch) < 1:
        parser.error("Use at least 100 pilot/evaluation replicas and positive counts/batch")
    if args.full_count % 2:
        parser.error("Use an even full count so the equal-budget signed alternative has exactly the same draw budget")
    args.output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).read_bytes()
    (args.output / "frozen_mixing_source.py").write_bytes(source)
    summary = {"model": "f=(1-epsilon) N(0,1)+epsilon N(2,1); M=N(0,1); g=epsilon(N(2,1)-N(0,1))",
        "observables": list(OBSERVABLES), "frequency": 1., "shift": 2.,
        "arguments": {k: str(v) if isinstance(v, Path) else v for k,v in vars(args).items()},
        "python": platform.python_version(), "numpy": np.__version__,
        "source_sha256": hashlib.sha256(source).hexdigest(), "cases": []}
    streams = np.random.SeedSequence(args.seed).spawn(4*9)
    for case, epsilon in enumerate((.05, .25, .5, .9)):
        print(f"epsilon={epsilon}", flush=True)
        seeds = streams[case*9:(case+1)*9]
        truth, baseline, vf, vd = moments(epsilon, 2., 1.)
        pilot_start = time.perf_counter()
        pf, pd, ptf, ptd = estimates(args.pilot, args.full_count, args.sign_count,
                                    epsilon, 2., 1., seeds[:3], args.batch)
        pvf, pvd = pf.var(axis=0, ddof=1), pd.var(axis=0, ddof=1)
        pc = np.mean((pf-pf.mean(axis=0))*(pd-pd.mean(axis=0)), axis=0)*args.pilot/(args.pilot-1)
        weight = np.clip((pvd-pc)/(pvf+pvd-2*pc), 0, 1)
        pilot_time = time.perf_counter()-pilot_start
        full, signed, tf, td = estimates(args.repetitions, args.full_count, args.sign_count,
                                       epsilon, 2., 1., seeds[3:6], args.batch)
        start = time.perf_counter()
        mixed = weight*full+(1-weight)*signed
        tm = time.perf_counter()-start
        # Both single-method alternatives spend the mixture's total Gaussian-draw budget.
        budget = args.full_count + 2*args.sign_count
        bf, bd, tbf, tbd = estimates(args.repetitions, budget, budget//2,
                                    epsilon, 2., 1., seeds[6:9], args.batch)
        raw = args.output / f"replicas_epsilon_{epsilon:g}.csv"
        with raw.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow([f"{method}_{obs}" for method in ("full", "signed", "mixed", "budget_full", "budget_signed") for obs in OBSERVABLES])
            writer.writerows(np.column_stack((full, signed, mixed, bf, bd)))
        for j, observable in enumerate(OBSERVABLES):
            vfull, vsigned = vf[j]/args.full_count, vd[j]/args.sign_count
            ideal = vsigned/(vfull+vsigned)
            covariance = float(np.cov(full[:,j], signed[:,j], ddof=1)[0,1])
            methods = {name: stats(values[:,j], truth[j]) for name,values in
                       (("full",full),("signed",signed),("mixed",mixed),("budget_full",bf),("budget_signed",bd))}
            summary["cases"].append({"epsilon": epsilon, "observable": observable, "truth": float(truth[j]),
                "pilot_full_weight": float(weight[j]), "exact_full_weight": float(ideal),
                "pilot_full_variance": float(pvf[j]), "pilot_signed_variance": float(pvd[j]), "pilot_covariance": float(pc[j]),
                "evaluation_covariance": covariance, "evaluation_correlation": float(np.corrcoef(full[:,j],signed[:,j])[0,1]),
                "exact_full_variance": float(vfull), "exact_signed_variance": float(vsigned),
                "exact_ideal_mixed_variance": float(vfull*vsigned/(vfull+vsigned)),
                "exact_pilot_weight_mixed_variance": float(weight[j]**2*vfull+(1-weight[j])**2*vsigned),
                "methods": methods, "runtime_seconds_all_observables": {"full": tf,"signed":td,"mix_arithmetic":tm,
                    "mixed_from_scratch":tf+td+tm,"pilot":pilot_time,"budget_full":tbf,"budget_signed":tbd},
                "both_available_variance_gain": min(methods['full']['variance'],methods['signed']['variance'])/methods['mixed']['variance'],
                "equal_draw_budget_variance_gain": min(methods['budget_full']['variance'],methods['budget_signed']['variance'])/methods['mixed']['variance'],
                "from_scratch_mse_time_gain": min(methods['budget_full']['mse']*tbf,methods['budget_signed']['mse']*tbd)/(methods['mixed']['mse']*(tf+td+tm)),
                "from_scratch_mse_time_gain_with_pilot": min(methods['budget_full']['mse']*tbf,methods['budget_signed']['mse']*tbd)/(methods['mixed']['mse']*(tf+td+tm+pilot_time))})
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    lines = ["# Frozen-distribution mixing experiment", "", summary["model"], "",
        f"{args.pilot} independent pilot replicas select observable-specific covariance-aware weights; {args.repetitions} fresh evaluation replicas measure performance.",
        "The Maxwellian observable is integrated analytically. Full and signed sample streams are independent. All expectations and variances are known exactly.", "",
        "Gain > 1 favors mixing. 'Both available' assumes component samples have already been generated. 'Equal draw budget' charges for generating both, and compares against each single estimator using that total Gaussian-draw budget.", "",
        "| epsilon | Observable | Pilot/exact full weight | Both-available variance gain | Equal-draw-budget variance gain | Measured/predicted mixed variance |", "|---:|---|---:|---:|---:|---:|"]
    for row in summary["cases"]:
        lines.append(f"| {row['epsilon']:g} | {row['observable']} | {row['pilot_full_weight']:.4f}/{row['exact_full_weight']:.4f} | {row['both_available_variance_gain']:.3f} | {row['equal_draw_budget_variance_gain']:.3f} | {row['methods']['mixed']['variance']/row['exact_pilot_weight_mixed_variance']:.3f} |")
    lines += ["", "## Bias and covariance", "", "| epsilon | Observable | Mixed bias / mean SE | Full-signed correlation |", "|---:|---|---:|---:|"]
    for row in summary["cases"]:
        m = row["methods"]["mixed"]
        lines.append(f"| {row['epsilon']:g} | {row['observable']} | {m['bias']/m['mean_standard_error']:.3f} | {row['evaluation_correlation']:.4f} |")
    lines += ["", "## Runtime", "", "Times include sampling and all three observable evaluations; they are not per-observable costs. Pilot cost is listed separately and must be charged unless amortized over repeated use.", "",
        "| epsilon | Full (s) | Signed (s) | Both + blend (s) | Blend arithmetic (s) | Equal-budget full/signed (s) | Pilot (s) |", "|---:|---:|---:|---:|---:|---:|---:|"]
    for row in summary["cases"][::3]:
        t = row["runtime_seconds_all_observables"]
        lines.append(f"| {row['epsilon']:g} | {t['full']:.3f} | {t['signed']:.3f} | {t['mixed_from_scratch']:.3f} | {t['mix_arithmetic']:.6f} | {t['budget_full']:.3f}/{t['budget_signed']:.3f} | {t['pilot']:.3f} |")
    lines += ["", "## Measured efficiency from scratch", "",
        "Gain is best single-estimator MSE times runtime divided by mixed MSE times runtime. Gain > 1 favors mixing. This charges generation of both components and includes all three observable evaluations in each runtime. Measurements are preliminary single timing blocks, not solver timings.", "",
        "| epsilon | Observable | MSE-time gain, pilot amortized | MSE-time gain, pilot charged |", "|---:|---|---:|---:|"]
    for row in summary['cases']:
        lines.append(f"| {row['epsilon']:g} | {row['observable']} | {row['from_scratch_mse_time_gain']:.3f} | {row['from_scratch_mse_time_gain_with_pilot']:.3f} |")
    lines += ["", "## Interpretation", "",
        "This controlled test isolates unbiased scalar mixing. It does not test the C++ population-count proxy, evolving fields, source sampling, projection, conservation, or shared-particle covariance.",
        "At fixed total sampling budget, independent estimators with variance proportional to inverse count favor allocating all draws to the estimator with lower variance per cost; mixing alone cannot beat that ideal allocation. Existing dual representations can still benefit from blending at negligible incremental cost.",
        "Raw replicas, exact predictions, measured MSE, covariance, weights, runtime and source provenance are saved. Timing is vectorized NumPy performance on this machine, not a C++ solver speedup."]
    (args.output / "REPORT.md").write_text("\n".join(lines)+"\n")
    # A statistical consistency check, with generous simultaneous tolerances.
    for row in summary["cases"]:
        ratio = row["methods"]["mixed"]["variance"]/row["exact_pilot_weight_mixed_variance"]
        tolerance = max(.1, 8*np.sqrt(2/(args.repetitions-1)))
        if not 1-tolerance < ratio < 1+tolerance:
            raise RuntimeError(f"Variance consistency failed: {row['epsilon']} {row['observable']}: {ratio}")
    print(args.output / "REPORT.md")


if __name__ == "__main__":
    main()
