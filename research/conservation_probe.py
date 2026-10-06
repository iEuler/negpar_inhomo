"""Test the existing signed-moment conservation option on adaptive pilot cases."""
import argparse
import json
from pathlib import Path
from benchmark import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot", type=Path)
    parser.add_argument("refinement", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    executable = (args.pilot / "provenance/negpar_inhomo.exe").resolve()
    results = []
    for label in ("base", "space"):
        source = args.refinement / f"adaptive_{label}"
        settings = json.loads((source / "input.json").read_text())
        settings["resampling"] = {"conserve_weighted_moments": True}
        print(label, flush=True)
        try:
            result = run(executable, args.output / label, settings, 30000, 600)
            result.update(status="complete", case=label)
        except Exception as error:
            result = {"status": "failed", "case": label, "failure": str(error)}
        results.append(result)
        (args.output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    lines = ["# Signed-moment conservation probe", "",
        "One seed, same frozen executable and settings as the refinement baseline,",
        "with conserve_weighted_moments enabled. This is an experimental option.", "",
        "| Grid | Status | Maximum relative energy drift |", "|---|---|---:|"]
    for result in results:
        value = f"{result['maximum_relative_energy_drift']:.4g}" if result["status"] == "complete" else "unavailable"
        lines.append(f"| {result['case']} | {result['status']} | {value} |")
    (args.output / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
