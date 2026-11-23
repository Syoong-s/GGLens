#!/usr/bin/env python3
"""Compare tangential ESD measurements relative to the first input table."""

import argparse
import warnings

import numpy as np

from esd_plot import configure_ratio_axis, display_radius


# ==========================================
# Function: parse_args
# Method: Match ordered input/label lists and validate axis parameters.
# ==========================================
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", required=True, nargs="+", action="extend",
        help="Existing ESD result tables; the first input is the reference. "
             "May be repeated or listed together; rows are paired by position.",
    )
    parser.add_argument(
        "--label", required=True, nargs="+", action="extend",
        help="One legend label per input, in the same order; quote spaces.",
    )
    parser.add_argument("--output", required=True, help="Output image path.")
    parser.add_argument("--title", default="BCG-centered group lensing")
    parser.add_argument(
        "--x-spread", type=float, default=0.06,
        help="Total display-only horizontal spread in log10(R) (default: 0.06).",
    )
    parser.add_argument(
        "--linear-width", type=float, default=0.1,
        help="Scale of the finely resolved region around ratio=1 (default: 0.1).",
    )
    parser.add_argument(
        "--ylim", type=float, nargs=2, metavar=("LOW", "HIGH"),
        help="Optional limits expressed as the original ratio values.",
    )
    args = parser.parse_args()
    if len(args.input) != len(args.label):
        parser.error("--input and --label must contain the same number of entries")
    if any(not label.strip() for label in args.label):
        parser.error("each input requires a nonempty label")
    if not np.isfinite(args.x_spread) or args.x_spread < 0:
        parser.error("--x-spread must be finite and nonnegative")
    if not np.isfinite(args.linear_width) or args.linear_width <= 0:
        parser.error("--linear-width must be finite and positive")
    if args.ylim is not None:
        if not np.all(np.isfinite(args.ylim)) or args.ylim[0] >= args.ylim[1]:
            parser.error("--ylim requires two finite values with LOW < HIGH")
    return args


# ==========================================
# Function: load_esd
# Method: Validate current and historical ESD tables while preserving every original row.
# ==========================================
def load_esd(path):
    table = np.loadtxt(path, comments="#", ndmin=2)
    if table.size == 0 or table.shape[1] not in (8, 11):
        raise ValueError(f"{path}: expected a nonempty 8- or 11-column ESD table")
    return table


# ==========================================
# Function: relative_ratio
# Method: Divide each original row by the same row of the fixed first input.
# ==========================================
def relative_ratio(table, reference, path, reference_path):
    if table.shape[0] != reference.shape[0]:
        raise ValueError(f"{path}: row count does not match first input {reference_path}")
    radius, measured, sigma, baseline = (
        table[:, 2], table[:, 4], table[:, 5], reference[:, 4]
    )
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        ratio = measured / baseline
    valid = (
        np.isfinite(radius) & (radius > 0) & np.isfinite(measured)
        & np.isfinite(baseline) & (baseline != 0) & np.isfinite(ratio)
    )
    if not np.any(valid):
        raise ValueError(f"{path}: no valid tangential ratio points")
    if not np.all(valid):
        warnings.warn(f"{path}: skipped {np.count_nonzero(~valid)} invalid rows")
    radius, ratio = radius[valid], ratio[valid]
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        # Treat the shared reference as fixed; show each input's own sigma.
        error = sigma[valid] / np.abs(baseline[valid])
        lower, upper = ratio - error, ratio + error
    valid_error = (
        np.isfinite(error) & (error >= 0)
        & np.isfinite(lower) & np.isfinite(upper)
    )
    if not np.all(valid_error):
        warnings.warn(
            f"{path}: {np.count_nonzero(~valid_error)} points have invalid "
            "uncertainties; keeping their ratios without error bars"
        )
    error = np.where(valid_error, error, 0)
    return radius, ratio, error


# ==========================================
    # Function: format_ratio
    # Method: Label the original ratio, never the transformed coordinate.
    # ==========================================
    def format_ratio(value, _position):
        return f"{value:.6g}"

    ax.yaxis.set_major_formatter(FuncFormatter(format_ratio))


# ==========================================
# Function: main
# Method: Normalize by the first input and emphasize even-numbered inputs.
# ==========================================
def main():
    args = parse_args()
    tables = [load_esd(path) for path in args.input]
    datasets = [
        relative_ratio(table, tables[0], path, args.input[0])
        for path, table in zip(args.input, tables)
    ]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))
    try:
        colors = plt.get_cmap("tab20")
        markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "*")
        for index, (label, data) in enumerate(zip(args.label, datasets)):
            radius, ratio, error = data
            emphasized = (index + 1) % 2 == 0
            ax.errorbar(
                display_radius(radius, index, len(datasets), args.x_spread),
                ratio, yerr=error,
                fmt=markers[index % len(markers)], color=colors((index ^ 1) % 20),
                markersize=6 if emphasized else 5, capsize=3,
                alpha=1.0 if emphasized else 0.5,
                zorder=3 if emphasized else 2, label=label,
            )
        ax.axhline(1, color="k", linestyle="--", linewidth=0.8)
        ax.set_xscale("log")
        configure_ratio_axis(ax, args.linear_width, args.ylim)
        ax.set_xlabel(r"$R\ [h^{-1}\,\mathrm{Mpc}]$")
        ax.set_ylabel(f"Measured / {args.label[0]} (first input)")
        ax.set_title("Tangential ratio")
        ax.grid(True, which="major", alpha=0.25)
        ax.legend()
        fig.suptitle(args.title)
        fig.tight_layout()
        fig.savefig(args.output, dpi=250)
    finally:
        plt.close(fig)


if __name__ == "__main__":
    main()
