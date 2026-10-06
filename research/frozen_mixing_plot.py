"""Render the controlled frozen mixing measurements with exact predictions."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    args = parser.parse_args()
    data = json.loads(args.summary.read_text())
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    labels = {"velocity": "Velocity", "kinetic_second_moment": "Second moment", "cosine_mode": "Cosine mode"}
    for observable, label in labels.items():
        rows = [r for r in data["cases"] if r["observable"] == observable]
        x = [r["epsilon"] for r in rows]
        line, = axes[0].plot(x, [r["both_available_variance_gain"] for r in rows], "o-", label=label)
        exact = [min(r["exact_full_variance"], r["exact_signed_variance"])/r["exact_pilot_weight_mixed_variance"] for r in rows]
        axes[0].plot(x, exact, ":", color=line.get_color(), alpha=.8)
        axes[1].plot(x, [r["equal_draw_budget_variance_gain"] for r in rows], "o-", label=label)
    axes[0].set_title("Both component estimates already available")
    axes[1].set_title("Both components sampled from scratch")
    for axis in axes:
        axis.axhline(1, color="black", linewidth=1, linestyle="--")
        axis.set_xlabel("Mixture fraction epsilon")
        axis.set_ylabel("Variance gain relative to best single estimator")
        axis.grid(alpha=.2)
        axis.legend(loc="lower right")
    axes[0].text(.02, .96, "Dotted: exact variance prediction\nGain > 1 favors mixing", transform=axes[0].transAxes, va="top", fontsize=9)
    axes[1].text(.02, .95, "Equal total Gaussian-draw budget\nPilot cost excluded; timings reported separately", transform=axes[1].transAxes, va="top", fontsize=9)
    figure.suptitle("Frozen Gaussian mixture: independent pilot weights, fresh evaluation samples")
    output = args.summary.parent / "variance_gain.png"
    figure.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
