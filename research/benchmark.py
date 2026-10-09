"""Reproducible Landau-damping pilot; never infer matched accuracy from weights.

Uses only the Python standard library. Raw runs, configurations, logs, wall
times and checks are retained. The reference is an ensemble of finer PIC runs,
not an exact solution. Use --seeds and --reference-weight for convergence work.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
import statistics
import subprocess
import shutil
import zipfile
import time


def series(directory, name):
    values = [float(v) for v in (directory / (name + ".txt")).read_text().split()]
    if not values or not all(math.isfinite(v) for v in values):
        raise ValueError(f"Non-finite or empty output: {directory / name}")
    return values


def config(mode, alpha, collision, weight, args):
    return {
        "schema_version": 1,
        "simulation": {
            "method": "pic" if mode == "pic" else "hdp",
            "spatial_cells": args.cells, "domain_length": 2 * math.pi,
            "time_step": args.dt, "signed_particle_weight": weight,
            "full_particle_weight": weight, "landau_amplitude": alpha,
            "collision_coefficient": collision, "poisson_coefficient": 1,
            "fourier_modes": args.modes,
        },
        "features": {
            "weighted_hdp": mode in ("mixed", "adaptive"),
            "weighted_fourier_resampling": mode in ("mixed", "adaptive"),
            "adaptive_effective_weights": mode == "adaptive",
        },
        "runtime": {"steps": args.steps, "threads": 1},
    }


def run(executable, output, settings, seed, timeout):
    output.mkdir(parents=True, exist_ok=False)
    settings["runtime"].update(seed=seed, output_directory=str(output.resolve()))
    cfg = output / "input.json"
    cfg.write_text(json.dumps(settings, indent=2) + "\n")
    start = time.perf_counter()
    with (output / "console.log").open("w") as log:
        completed = subprocess.run([str(executable), "--config", str(cfg)],
                                   stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
    wall = time.perf_counter() - start
    if completed.returncode:
        raise RuntimeError(f"Run failed ({completed.returncode}); see {output / 'console.log'}")
    # Validate every numerical output, including snapshots and distributions.
    for path in output.glob("*.txt"):
        if path.name == "run_metadata.txt":
            continue
        for token in path.read_text().split():
            if re.match(r"(?i)^[+-]?(?:nan(?:\([^)]*\))?|inf(?:inity)?|1\.#(?:inf|qnan|snan|ind))$", token):
                raise ValueError(f"Non-finite output in {path}")
            try:
                number = float(token)
            except ValueError:
                continue
            if not math.isfinite(number):
                raise ValueError(f"Non-finite output in {path}")
    name = "elec_energy_F" if settings["simulation"]["method"] == "pic" else "elec_energy"
    electric = series(output, name)
    times = series(output, "time_rec")
    expected = settings["runtime"]["steps"] + 1
    if len(electric) != expected or len(times) != expected:
        raise ValueError(f"Incomplete state history in {output}")
    record = {"directory": str(output), "seed": seed, "wall_seconds": wall,
              "electric_l2_squared": electric, "times": times,
              "cpu_seconds": sum(series(output, "cputime_all")),
              "final_signed_count": series(output, "Np_rec")[-1] + series(output, "Nn_rec")[-1],
              "final_full_count": series(output, "Nf_rec")[-1],
              "final_signed_weight": series(output, "Neff_D_rec")[-1],
              "final_full_weight": series(output, "Neff_F_rec")[-1]}
    energy_name = "total_energy_F" if settings["simulation"]["method"] == "pic" else "totalEnergy"
    energy = series(output, energy_name)
    record["maximum_relative_energy_drift"] = max(abs(v / energy[0] - 1) for v in energy)
    for name in ("mass_rec", "mass_F_rec"):
        if (output / (name + ".txt")).exists():
            mass = series(output, name)
            record[name + "_maximum_relative_drift"] = max(abs(v / mass[0] - 1) for v in mass)
    (output / "measurement.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--cells", type=int, default=40)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--modes", type=int, default=8)
    parser.add_argument("--weights", type=float, nargs="+", default=[0.001, 0.0005])
    parser.add_argument("--reference-weight", type=float, default=0.0001)
    parser.add_argument("--alphas", type=float, nargs="+", default=[0.05, 0.4])
    parser.add_argument("--collisions", type=float, nargs="+", default=[1.0, 10.0])
    parser.add_argument("--timeout", type=float, default=600)
    args = parser.parse_args()
    args.executable = args.executable.resolve()
    if args.seeds < 2 or args.steps < 1:
        parser.error("Use at least two independent seeds and one step")
    if args.output.exists():
        parser.error("Output directory must be new; previous results are never overwritten")
    args.output.mkdir(parents=True)
    provenance = args.output / "provenance"
    provenance.mkdir()
    source_root = Path(__file__).resolve().parents[1]
    with zipfile.ZipFile(provenance / "source.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for folder in ("src", "config", "tests", "research"):
            for path in (source_root / folder).rglob("*"):
                if path.is_file() and "runs" not in path.parts and "__pycache__" not in path.parts:
                    archive.write(path, path.relative_to(source_root))
        for name in ("CMakeLists.txt", "CMakePresets.json", "README.md"):
            archive.write(source_root / name, name)
    frozen_executable = provenance / args.executable.name
    shutil.copy2(args.executable, frozen_executable)
    for dll in args.executable.parent.glob("*.dll"):
        shutil.copy2(dll, provenance / dll.name)
    args.executable = frozen_executable.resolve()
    records = []
    summary = {"status": "running", "arguments": {k: str(v) if isinstance(v, Path) else v
              for k, v in vars(args).items()},
              "executable_sha256": hashlib.sha256(args.executable.read_bytes()).hexdigest(),
              "observable": "integral E(x)^2 dx (no square root)",
              "reference": "independent fine PIC ensemble; not an exact solution", "cases": []}
    summary_path = args.output / "summary.json"
    try:
        for alpha in args.alphas:
            for collision in args.collisions:
                label = f"alpha_{alpha:g}_A_{collision:g}"
                references = []
                for seed in range(args.seeds):
                    print(f"{label} reference seed {seed}", flush=True)
                    references.append(run(args.executable, args.output / label / f"reference_{seed}",
                                          config("pic", alpha, collision, args.reference_weight, args),
                                          10000 + seed, args.timeout))
                mean_reference = [statistics.mean(row) for row in zip(
                    *(r["electric_l2_squared"] for r in references))]
                for mode in ("pic", "hdp", "mixed", "adaptive"):
                    for weight in args.weights:
                        runs = []
                        for seed in range(args.seeds):
                            print(f"{label} {mode} weight {weight:g} seed {seed}", flush=True)
                            runs.append(run(args.executable,
                                args.output / label / f"{mode}_w_{weight:g}_seed_{seed}",
                                config(mode, alpha, collision, weight, args), 20000 + seed, args.timeout))
                        normalizer = sum(v * v for v in mean_reference)
                        errors = [math.sqrt(sum((a-b)**2 for a, b in zip(
                            r["electric_l2_squared"], mean_reference)) / normalizer) for r in runs]
                        records.append({"case": label, "mode": mode, "weight": weight,
                            "trajectory_relative_rms_error_mean": statistics.mean(errors),
                            "trajectory_relative_rms_error_sd": statistics.stdev(errors),
                            "wall_seconds_mean": statistics.mean(r["wall_seconds"] for r in runs),
                            "wall_seconds_sd": statistics.stdev(r["wall_seconds"] for r in runs),
                            "runs": runs})
                        summary["cases"] = records
                        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
        summary["status"] = "complete"
    except Exception as error:
        summary["status"] = "failed"
        summary["failure"] = str(error)
        raise
    finally:
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
