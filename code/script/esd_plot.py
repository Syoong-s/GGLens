"""Plotting utilities for the fitted ESD and optional table comparison."""

from contextlib import ExitStack
import fcntl
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory

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
# Method: Draw the best-fit components and total posterior bands against the measured ESD.
# ==========================================
def render_fit(radius, measured, covariance, profiles, chi2, dof, path, title):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    r = profiles.radius
    best = profiles.best[:, -1]
    p025, p16, _, p84, p975 = profiles.quantiles[:, :, -1]
    fig, ax = plt.subplots(figsize=(7.4, 5.1))
    ax.axhline(0, color="0.62", linewidth=0.8, zorder=0)
    ax.fill_between(r, p025, p975, color="#c9ddf1", label="95% posterior")
    ax.fill_between(r, p16, p84, color="#83b0dc", label="68% posterior")
    ax.plot(r, best, color="#123b65", linewidth=2, label="Total best fit")
    component_styles = (
        ("Stellar", "#b6542b", ":"),
        ("Centered one-halo", "#247a69", "--"),
        ("Off-center one-halo", "#b14f61", (0, (5, 2, 1, 2, 1, 2))),
        ("Two-halo", "#80478a", "-."),
    )
    for values, (label, color, linestyle) in zip(profiles.best.T[:-1], component_styles):
        ax.plot(r, values, label=label, color=color, linestyle=linestyle, linewidth=1.4)
    ax.errorbar(radius, measured, yerr=np.sqrt(np.diag(covariance)),
                color="#252525", fmt="o", ms=4, capsize=2, zorder=3, label="Tangential ESD")
    ax.set(xscale="log", yscale="symlog", xlabel=r"$R\ [h^{-1}\,\mathrm{Mpc}]$",
           ylabel=r"$\Delta\Sigma\ [h\,M_\odot/\mathrm{pc}^2]$", title=title)
    ax.text(0.03, 0.06, rf"$\chi^2/\nu={chi2:.2f}/{dof}$", transform=ax.transAxes)
    ax.legend(frameon=False, ncol=2, fontsize=9)
    save_figure(fig, path, 220)


def render_fd_test(path, points, calibration, limit):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0))
    fig.suptitle("Field-distortion self-calibration")
    field = np.array([-limit, limit])
    for component, (ax, rows, (m, c)) in enumerate(zip(axes, points, calibration), 1):
        gf, shear, sigma = np.asarray(rows).reshape(-1, 3).T
        ax.plot(field, field, color="0.55", linestyle="--", linewidth=1.1, label=r"$\hat g = g_f$")
        ax.plot(field, (1 + m)*field + c, color="#247a69", linewidth=1.5,
                label=rf"Fit  $m={m:.3g}$, $c={c:.3g}$")
        ax.errorbar(gf, shear, yerr=sigma, fmt="o", color="#123b65", ms=3.5, capsize=2,
                    zorder=3, label="FD bins")
        ax.set(xlabel=rf"Field distortion $g_{{f,{component}}}$", title=rf"Component $g_{component}$")
        ax.ticklabel_format(style="sci", scilimits=(-2, 2), useMathText=True)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel(r"Measured shear $\hat g$")
    lower = min(ax.get_ylim()[0] for ax in axes)
    upper = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(lower, upper)
    save_figure(fig, path, 220)


def posterior_percentiles(samples, weights):
    order = np.argsort(samples)
    cumulative = np.cumsum(weights[order], dtype=np.int64)
    positions = (cumulative[-1] - 1) * np.array([.16, .5, .84])
    low = samples[order[np.searchsorted(cumulative, np.floor(positions), side="right")]]
    high = samples[order[np.searchsorted(cumulative, np.ceil(positions), side="right")]]
    fraction = positions - np.floor(positions)
    return np.where(fraction < .5, low + (high - low)*fraction,
                    high - (high - low)*(1 - fraction))


def render_posterior(samples, weights, path):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    labels = (r"$\log_{10} M_{180m}$", r"$c$", r"$f_{\rm cen}$", r"$R_{\rm sig}/r_{180m}$", r"$b_h$")
    size = len(labels)
    edges = [np.histogram_bin_edges(values, bins=36) for values in samples]
    fig, axes = plt.subplots(size, size, figsize=(10.5, 10))
    for i, label in enumerate(labels):
        for j in range(size):
            ax = axes[i, j]
            if j > i:
                ax.set_axis_off()
                continue
            if i == j:
                counts, _ = np.histogram(samples[i], bins=edges[i], weights=weights)
                ax.stairs(counts / counts.sum(), edges[i], fill=True, color="#83b0dc")
                low, median, high = posterior_percentiles(samples[i], weights)
                for location, style in ((low, "--"), (median, "-"), (high, "--")):
                    ax.axvline(location, color="#123b65", linestyle=style, linewidth=1)
                ax.set_title(f"{label}  {median:.4g} [{low:.4g}, {high:.4g}]", fontsize=9)
                ax.set_yticks([])
            else:
                density, _, _ = np.histogram2d(samples[j], samples[i], bins=(edges[j], edges[i]), weights=weights)
                ax.pcolormesh(edges[j], edges[i], density.T, cmap="Blues", shading="auto")
                if j == 0:
                    ax.set_ylabel(label)
                else:
                    ax.tick_params(labelleft=False)
            if i == size - 1:
                ax.set_xlabel(labels[j])
                ax.tick_params(axis="x", labelrotation=35)
            else:
                ax.tick_params(labelbottom=False)
    save_figure(fig, path, 180)


def read_legacy_fit(prefix, posterior=False):
    from esd_archive import PARAMETERS
    from esd_profiles import read_profiles, split_one_halo

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
        best = [summary["parameters"][name]["best"] for name in PARAMETERS]
        profiles = split_one_halo(profiles, best, summary["cosmology"]["omega_m"])
        samples = chain["posterior"].reshape(-1, 5).T if posterior else None
        return summary, profiles, (obs, chain["measured"], chain["covariance"]), samples


def render_saved_fit(prefix, path):
    summary, profiles, observed, _ = read_legacy_fit(prefix)
    render_fit(*observed, profiles, summary["chi2"], summary["degrees_of_freedom"], path, summary["plot_title"])


def open_saved_fit(archive, legacy_base):
    from esd_archive import load_fit

    archive = Path(archive)
    if archive.is_symlink() or archive.exists():
        fit = load_fit(archive)
        chain, _ = fit
        return fit, chain.measurement_path, chain.measurement_id
    if Path(f"{legacy_base}_fit.json").is_file():
        return read_legacy_fit(legacy_base, posterior=True), "", None
    raise FileNotFoundError(f"No completed fit archive at {archive}")


def publish_saved_fit(saved, paths, radii, measured, covariance, context, table):
    from esd_archive import lock_directory
    from esd_profiles import plot_profiles, split_one_halo

    radii = np.asarray(radii)
    measured = np.asarray(measured)
    covariance = np.asarray(covariance).reshape((len(radii), len(radii)))
    title = context["title"]
    if len(saved) == 2:
        chain, draws = saved
        _, best = chain.best_sample()
        chi2, dof = -2 * best[-1], len(radii) - 5
        profiles = plot_profiles(draws, radii, context["stellar_mass"])
        omega_m = chain.legacy_cosmology[0] if chain.legacy_cosmology.size else context["omega_m"]
        profiles = split_one_halo(profiles, best, omega_m)
        samples, weights = chain.states[:5], chain.weights
        title = title or chain.title or f'{context["center"].upper()}-centered group lensing'
    else:
        summary, profiles, observed, samples = saved
        if (context["center"] != summary["center"] or
                any(not np.array_equal(a, b) for a, b in zip(observed, (radii, measured, covariance)))):
            raise ValueError("Legacy fit and measurement contain different ESD data")
        chi2, dof = summary["chi2"], summary["degrees_of_freedom"]
        title = title or summary["plot_title"]
        weights = np.ones(samples.shape[1], dtype=np.int64)

    destinations = {name: Path(path) for name, path in paths.items()}
    with ExitStack() as cleanup:
        staged = {}
        for name, path in destinations.items():
            if path.is_dir():
                raise IsADirectoryError(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            directory = cleanup.enter_context(TemporaryDirectory(dir=path.parent, prefix=f".{path.name}."))
            staged[name] = Path(directory) / path.name
        staged["table"].write_bytes(table)
        render_fit(radii, measured, covariance, profiles, chi2, dof, staged["esd"], title)
        render_posterior(samples, weights, staged["posterior"])
        if "fd" in staged:
            render_fd_test(staged["fd"], **context["fd"])
        with lock_directory(destinations["table"].parent, fcntl.LOCK_EX):
            for name, path in destinations.items():
                os.replace(staged[name], path)
