"""Summarize completed pilot measurements without declaring unmatched speedups."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    args = parser.parse_args()
    data = json.loads(args.summary.read_text())
    if data["status"] != "complete":
        raise ValueError("Only complete experiments can be reported")
    output = args.summary.parent
    rows = data["cases"]
    header = ["# Reproducible Landau damping pilot", "",
        "These are pilot measurements against an independent fine-PIC ensemble, not exact errors.",
        "Error is relative RMS over the full trajectory of the squared electric-field norm.",
        "Runtime includes initialization, reconstruction, synchronization and output.",
        "No matched-error speedup is inferred from equal input particle weights.", "",
        f"Executable SHA-256: `{data['executable_sha256']}`", "",
        "| Case | Method | Input weight | Relative error, mean ± SD | Wall seconds, mean ± SD | Final signed/full count, mean |",
        "|---|---|---:|---:|---:|---:|"]
    latex = [r"\begin{tabular}{llrrr}", r"\hline",
        r"$(\alpha,A)$ & Method & Weight & Rel. RMS error & Wall time (s) \\", r"\hline"]
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["case"]].append(row)
        header.append(f"| {row['case']} | {row['mode']} | {row['weight']:g} | "
            f"{row['trajectory_relative_rms_error_mean']:.4g} ± {row['trajectory_relative_rms_error_sd']:.2g} | "
            f"{row['wall_seconds_mean']:.3g} ± {row['wall_seconds_sd']:.2g} | "
            f"{statistics.mean(r['final_signed_count'] for r in row['runs']):.0f}/"
            f"{statistics.mean(r['final_full_count'] for r in row['runs']):.0f} |")
        alpha = row["case"].split("_")[1]
        collision = row["case"].split("_")[3]
        latex.append(f"$({alpha},{collision})$ & {row['mode']} & {row['weight']:g} & "
            f"{row['trajectory_relative_rms_error_mean']:.3g} & {row['wall_seconds_mean']:.3g} " + r"\\")
    header += ["", "Interpretation limits:", "",
        "- Three seeds provide a preliminary spread estimate, not a precise confidence interval.",
        "- The reference uses the same grid and timestep; it cannot detect shared discretization bias.",
        "- The plotted cases start with a local Maxwellian. They do not cover arbitrary far-from-equilibrium initial data.",
        "- A Fourier cutoff of eight is a pilot setting and requires refinement.",
        "- The source sampler is the currently implemented legacy construction.",
        "- Adaptive constraints preserve a population-count variance proxy; full observable errors are measured independently.",
        "- Conservation and reference-refinement checks remain necessary before publication."]
    header += ["", "| Case | Method | Weight | Maximum relative energy drift, mean |", "|---|---|---:|---:|"]
    for row in rows:
        header.append(f"| {row['case']} | {row['mode']} | {row['weight']:g} | "
                      f"{statistics.mean(r['maximum_relative_energy_drift'] for r in row['runs']):.4g} |")
    (output / "REPORT.md").write_text("\n".join(header) + "\n", encoding="utf-8")
    latex += [r"\hline", r"\end{tabular}"]
    (output / "pilot_table.tex").write_text("\n".join(latex) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    colors = {"pic": "#444444", "hdp": "#2474b5", "mixed": "#dd8523", "adaptive": "#19814a"}
    for axis, (label, values) in zip(axes.flat, sorted(grouped.items())):
        for mode in colors:
            group = sorted((r for r in values if r["mode"] == mode), key=lambda r: r["weight"])
            axis.errorbar([r["trajectory_relative_rms_error_mean"] for r in group],
                [r["wall_seconds_mean"] for r in group],
                xerr=[r["trajectory_relative_rms_error_sd"] for r in group],
                yerr=[r["wall_seconds_sd"] for r in group],
                fmt="o-", capsize=3, color=colors[mode], label=mode.upper())
        axis.set_title(label.replace("_", " "))
        axis.set_xlabel("Relative trajectory RMS error (fine-PIC ensemble)")
        axis.set_ylabel("Total wall time (s)")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, alpha=0.2)
        axis.legend(fontsize=8)
    figure.suptitle("Pilot cost versus achieved error; three independent seeds", fontsize=13)
    figure.savefig(output / "cost_error.png", dpi=180)
    print(output / "REPORT.md")
    print(output / "cost_error.png")


if __name__ == "__main__":
    main()
