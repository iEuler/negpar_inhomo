"""Small seeded refinement and exact replay checks using a frozen pilot solver."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
from benchmark import config, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    executable = (args.pilot / "provenance/negpar_inhomo.exe").resolve()
    args.output.mkdir(parents=True, exist_ok=args.resume)
    records = []
    # Equal final time; these are preliminary sensitivity checks, not convergence orders.
    for mode in ("pic", "hdp", "adaptive"):
        for label, cells, dt, steps, modes, weight in (
                ("base", 40, .01, 100, 8, .0005),
                ("time", 40, .005, 200, 8, .0005),
                ("space", 80, .01, 100, 8, .0005),
                ("fourier", 40, .01, 100, 12, .0005),
                ("particles", 40, .01, 100, 8, .00025)):
            settings = config(mode, .4, 1, weight,
                SimpleNamespace(cells=cells, dt=dt, steps=steps, modes=modes))
            print(mode, label, flush=True)
            directory = args.output / f"{mode}_{label}"
            try:
                if args.resume and directory.exists():
                    result = json.loads((directory / "measurement.json").read_text())
                else:
                    result = run(executable, directory, settings, 30000, 600)
            except Exception as error:
                records.append({"mode": mode, "refinement": label, "status": "failed",
                    "failure": str(error), "directory": str(directory)})
                continue
            result["status"] = "complete"
            result.update(mode=mode, refinement=label)
            records.append(result)
        baseline = args.output / f"{mode}_base"
        replay = args.output / f"{mode}_replay"
        settings = json.loads((baseline / "input.json").read_text())
        if not (args.resume and (replay / "measurement.json").exists()):
            run(executable, replay, settings, 30000, 600)
        files = [p for p in baseline.glob("*.txt")
                 if p.name != "run_metadata.txt" and not p.name.startswith("cputime")]
        mismatches = [p.name for p in files if p.read_bytes() != (replay / p.name).read_bytes()]
        if mismatches:
            raise RuntimeError(f"Seeded replay differs: {mode}: {mismatches}")
    lines = ["# Preliminary refinement and replay checks", "",
        "One seed, alpha=0.4, A=1, final time=1. Sensitivity is not a convergence proof.", "",
        f"Frozen executable SHA-256: `{hashlib.sha256(executable.read_bytes()).hexdigest()}`", "",
        "PIC, HDP and adaptive numerical outputs replay identically with the same seed.", "",
        "| Mode | Refinement | Relative trajectory change from baseline | Maximum relative energy drift |", "|---|---|---:|---:|"]
    for row in records:
        if row["status"] == "failed":
            lines.append(f"| {row['mode']} | {row['refinement']} | FAILED; see console.log | unavailable |")
            continue
        base = next(r for r in records if r["mode"] == row["mode"] and r["refinement"] == "base")
        values = row["electric_l2_squared"][::2] if row["refinement"] == "time" else row["electric_l2_squared"]
        reference = base["electric_l2_squared"]
        delta = math.sqrt(sum((a-b)**2 for a, b in zip(values, reference)) / sum(b*b for b in reference))
        row["trajectory_change_from_baseline"] = delta
        lines.append(f"| {row['mode']} | {row['refinement']} | {delta:.4g} | {row['maximum_relative_energy_drift']:.4g} |")
    (args.output / "summary.json").write_text(json.dumps(records, indent=2) + "\n")
    (args.output / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
