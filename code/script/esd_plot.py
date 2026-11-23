"""Plotting utilities for the fitted ESD and optional table comparison."""

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
# Function: render_fit
# Method: Draw the saved total best-fit curve and posterior bands against the measured ESD.
# ==========================================
def render_fit(radius, measured, covariance, profiles, chi2, dof, path, title):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    r = profiles.radius
    best = profiles.best[:, -1]
    p025, p16, _, p84, p975 = profiles.quantiles[:, :, -1]
    fig, ax = plt.subplots(figsize=(7.4, 5.1))
    ax.fill_between(r, p025, p975, color="#c9ddf1", label="95% posterior")
    ax.fill_between(r, p16, p84, color="#83b0dc", label="68% posterior")
    ax.plot(r, best, color="#123b65", linewidth=2, label="Best-fit model")
    ax.errorbar(radius, measured, yerr=np.sqrt(np.diag(covariance)),
                color="#252525", fmt="o", ms=4, capsize=2, zorder=3, label="Tangential ESD")
    ax.axhline(0, color="0.62", linewidth=0.8)
    ax.set(xscale="log", yscale="symlog", xlabel=r"$R\ [h^{-1}\,\mathrm{Mpc}]$",
           ylabel=r"$\Delta\Sigma\ [h\,M_\odot/\mathrm{pc}^2]$", title=title)
    ax.text(0.03, 0.06, rf"$\chi^2/\nu={chi2:.2f}/{dof}$", transform=ax.transAxes)
    ax.legend(frameon=False)
    save_figure(fig, path, 220)


# ==========================================
# Function: render_saved_fit
# Method: Reproduce the fit plot from one result prefix and its saved plot data.
# ==========================================
def render_saved_fit(prefix, path):
    import json
    from esd_profiles import read_profiles

    summary = json.loads(Path(f"{prefix}_fit.json").read_text())
    fit_id = summary["fit_id"]
    profiles = read_profiles(f"{prefix}_profiles.dat", fit_id)
    with np.load(f"{prefix}_chain.npz", allow_pickle=False) as chain:
        if "fit_id" not in chain or chain["fit_id"].item() != fit_id:
            raise ValueError("Fit summary and chain have different fit IDs")
        if chain["center"].item() != summary["center"]:
            raise ValueError("Fit summary and chain have different lens centers")
        obs = chain["radii"]
        if not np.array_equal(obs, np.asarray(summary["fitted_profile"]["R_Mpc_over_h"])):
            raise ValueError("Fit summary and chain have different measured radii")
        keep = np.searchsorted(profiles.radius, obs)
        if (np.any(keep >= len(profiles.radius)) or not np.array_equal(profiles.radius[keep], obs)
                or not np.array_equal(profiles.best[keep, -1], summary["fitted_profile"]["total"])):
            raise ValueError("Fit summary and plotted best-fit values disagree")
        render_fit(obs, chain["measured"], chain["covariance"], profiles, summary["chi2"],
                   summary["degrees_of_freedom"], path, summary["plot_title"])
