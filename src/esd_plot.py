"""Render ESD and Cartesian FD results passed in memory by the C++ program."""

import os
from pathlib import Path

import numpy as np


# ==========================================
# Function: display_radius
# Method: Center a bounded group spread in log10(R), without changing data.
# ==========================================
def display_radius(radius, index, count, spread):
    offset = 0.0 if count == 1 else spread * (index / (count - 1) - 0.5)
    return np.asarray(radius) * 10.0 ** offset


# ==========================================
# Function: configure_ratio_axis
# Method: Compress deviations with asinh and place ticks in original units.
# ==========================================
def configure_ratio_axis(ax, width, limits):
    # ==========================================
    # Function: forward
    # Method: Smooth linear central region and logarithmic outer compression.
    # ==========================================
    def forward(values):
        return np.arcsinh((np.asarray(values) - 1.0) / width)

    # ==========================================
    # Function: inverse
    # Method: Recover the original ratio values for Matplotlib transforms.
    # ==========================================
    def inverse(values):
        with np.errstate(over="ignore", invalid="ignore"):
            return 1.0 + width * np.sinh(values)

    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    ax.set_yscale("function", functions=(forward, inverse))
    if limits is not None:
        ax.set_ylim(*limits)
    low, high = ax.get_ylim()
    transformed_low, transformed_high = forward([low, high])
    min_spacing = (transformed_high - transformed_low) / 12.0
    extent = max(abs(low - 1.0), abs(high - 1.0))
    # Prefer ratio=1, then nearby rounded offsets, followed by coarser decades.
    offsets = [0.0]
    if extent > 0:
        first = int(np.floor(np.log10(width))) - 2
        last = int(np.ceil(np.log10(extent))) + 1
        offsets.extend(
            mantissa * 10.0 ** exponent
            for exponent in range(first, last + 1)
            for mantissa in (1.0, 2.0, 5.0)
        )
    ticks, positions = [], []
    for offset in offsets:
        candidates = (1.0,) if offset == 0 else (1.0 - offset, 1.0 + offset)
        for value in candidates:
            if not low <= value <= high:
                continue
            position = float(forward(value))
            if all(abs(position - old) >= min_spacing for old in positions):
                ticks.append(value)
                positions.append(position)
    if len(ticks) < 2:
        # Narrow manual limits may contain no rounded candidates.
        ticks = inverse(np.linspace(transformed_low, transformed_high, 5))
    ax.yaxis.set_major_locator(FixedLocator(sorted(ticks)))
    ax.yaxis.set_minor_locator(NullLocator())

    # ==========================================
# Function: save_figure
# Method: Commit a complete PNG in the destination directory and release figure state.
# ==========================================
def save_figure(fig, output, dpi):
    from matplotlib import pyplot as plt

    path = Path(output)
    temp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(temp, format="png", dpi=dpi)
        os.replace(temp, path)
    finally:
        plt.close(fig)
        temp.unlink(missing_ok=True)


# ==========================================
# Function: render_esd
# Method: Plot measured and NFW ESD, the cross component, and the tangential ratio.
# ==========================================
def render_esd(rows, output, title, x_spread, linear_width, limits):
    from matplotlib import pyplot as plt

    radius, measured, sigma, cross, cross_sigma, theory, ratio = np.asarray(rows).T
    valid = np.isfinite(radius) & (radius > 0) & np.isfinite(measured)
    valid_theory = valid & np.isfinite(theory)
    valid_cross = valid & np.isfinite(cross)
    ratio_valid = valid_theory & np.isfinite(ratio) & (theory != 0)
    # Shift the plotted components symmetrically on the logarithmic R axis.
    # Keep theory, ratio-panel radii and all measured values unmodified.
    tangential_radius = display_radius(radius, 0, 2, x_spread)
    cross_radius = display_radius(radius, 1, 2, x_spread)
    cross_order = np.flatnonzero(valid_cross)
    cross_order = cross_order[np.argsort(radius[cross_order], kind="stable")]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    ax.errorbar(tangential_radius[valid], measured[valid],
                yerr=np.where(np.isfinite(sigma[valid]) & (sigma[valid] >= 0), sigma[valid], 0),
                fmt="o", color="tab:red", capsize=3, zorder=3,
                label="Tangential PDF_SYM")
    ax.errorbar(cross_radius[cross_order], cross[cross_order],
                yerr=np.where(np.isfinite(cross_sigma[cross_order]) & (cross_sigma[cross_order] >= 0),
                              cross_sigma[cross_order], 0),
                fmt="s--", color="tab:blue", capsize=2,
                alpha=0.4, linewidth=1, zorder=1, label="Cross component")
    ax.plot(radius[valid_theory], theory[valid_theory], "k--", label="NFW reference")
    ax.axhline(0, color="0.5", linewidth=0.7)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel(r"$R\ [h^{-1}\,\mathrm{Mpc}]$")
    ax.set_ylabel(r"$\Delta\Sigma\ [h\,M_\odot/\mathrm{pc}^2]$")
    ax.set_title("ESD and cross component")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[1]
    ax.errorbar(radius[ratio_valid], ratio[ratio_valid],
                yerr=np.where(np.isfinite(sigma[ratio_valid]) & (sigma[ratio_valid] >= 0),
                              sigma[ratio_valid] / np.abs(theory[ratio_valid]), 0),
                fmt="o", color="tab:red", capsize=3)
    ax.axhline(1, color="k", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    configure_ratio_axis(ax, linear_width, limits)
    ax.set_xlabel(r"$R\ [h^{-1}\,\mathrm{Mpc}]$")
    ax.set_ylabel("Measured / NFW")
    ax.set_title("Tangential ratio")
    ax.grid(True, alpha=0.25)

    fig.suptitle(title)
    save_figure(fig, output, 250)


# ==========================================
# Function: render_fd
# Method: Plot both Cartesian components against FD with fitted gains and offsets.
# ==========================================
def render_fd(output, selected, limit, components):
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x = np.linspace(-limit, limit, 200)
    for component, (ax, fit) in enumerate(zip(axes, components), 1):
        m, c, sigma_m, sigma_c, n_fit, rows = fit
        points = np.asarray(rows, dtype=float)
        valid = ((points[:, 0] > 0) & (points[:, 4] == 0)
                 & np.all(np.isfinite(points), axis=1) & (points[:, 3] > 0))
        ax.errorbar(points[valid, 1], points[valid, 2], yerr=points[valid, 3],
                    fmt="o", markersize=4, capsize=2, label="All-R unique sources")
        ax.plot(x, x, ":", color="black", label="Ideal")
        ax.plot(x, (1 + m) * x + c, color="tab:red", label="Weighted fit")
        ax.set_xlabel(rf"$g_{{f,{component}}}$")
        ax.set_ylabel(rf"$\hat g_{{{component}}}$")
        ax.set_title(f"Cartesian component {component}")
        text = (f"m{component} = {m:.5g} ± {sigma_m:.3g}\n"
                f"c{component} = {c:.5g} ± {sigma_c:.3g}\n"
                f"Valid FD bins: {n_fit}/{len(points)}")
        ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top", fontsize=9)
        ax.grid(alpha=0.25)
        ax.legend(loc="lower right", fontsize=8)
    fig.suptitle(f"All-R source self-calibration — {selected:,} unique measurement rows")
    save_figure(fig, output, 160)


# ==========================================
# Function: render_plots
# Method: Use the noninteractive backend for the single end-of-run plotting call.
# ==========================================
def render_plots(rows, output, title, x_spread, linear_width, limits, fd):
    import matplotlib

    matplotlib.use("Agg")
    render_esd(rows, output, title, x_spread, linear_width, limits)
    if fd is not None:
        render_fd(*fd)
