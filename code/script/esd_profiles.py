"""Generate and restore the joint ESD curve draws."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from esd_archive import CurveArchive
from esd_nfw import RHO_CRITICAL, halo_radius, nfw_esd, nfw_scales

COMPONENTS = ("stellar", "one_halo", "two_halo", "total")
PERCENTILES = (2.5, 16, 50, 84, 97.5)
COLUMNS = ("best", "p2p5", "p16", "p50", "p84", "p97p5")


@dataclass(frozen=True)
class Profiles:
    radius: np.ndarray
    best: np.ndarray
    quantiles: np.ndarray


def stellar_profile(radius, stellar_mass):
    return stellar_mass / (np.pi*radius*radius*1e12)


def evaluate_profile_draws(problem, posterior, best_index, limit, seed):
    from esd_two_halo import TwoHalo

    count = min(limit, len(posterior))
    if count == len(posterior):
        indices = np.arange(count, dtype="<u4")
    else:
        rng = np.random.default_rng(np.random.SeedSequence([seed, 4796]))
        indices = np.sort(rng.choice(len(posterior), count, replace=False)).astype("<u4")

    smooth = np.geomspace(problem.r.min(), problem.r.max(), 180)
    radius = np.union1d(smooth, problem.r)
    extra_radii = np.setdiff1d(radius, problem.r).astype("<f8", copy=False)
    candidates = posterior[np.append(indices, best_index)]
    _, keep, mapping = np.unique(candidates.view(np.uint64), axis=0,
                                 return_index=True, return_inverse=True)
    mapping = mapping.reshape(-1)
    theta = candidates[keep]
    rows = mapping[:-1].astype("<u4")
    best_row = int(mapping[-1])

    model = problem.model
    two_halo = TwoHalo(radius, model.correlation, model.mean_density)
    curves = np.empty((2, len(radius), len(theta)), dtype="<f8")
    batch = 256
    for first in range(0, len(theta), batch):
        last = min(len(theta), first + batch)
        _, one, two = model.components(theta[first:last], radius, two_halo)
        curves[0, :, first:last] = one.T
        curves[1, :, first:last] = two.T
    if not np.all(np.isfinite(curves)):
        raise ValueError("Non-finite ESD profile in the retained posterior")
    return CurveArchive(len(radius), extra_radii, indices, rows, best_row, curves)


def curve_terms(curves, observed_radii, stellar_mass):
    radius = curves.radius(observed_radii)
    one, two = curves.components
    stellar = stellar_profile(radius, stellar_mass)[:, None]
    return radius, stellar, one, two, (stellar + one) + two


def restore_draws(curves, observed_radii, stellar_mass):
    """Restore the original (draw, radius, component) curves and global best curve."""
    radius, stellar, one, two, total = curve_terms(curves, observed_radii, stellar_mass)
    terms = (np.broadcast_to(stellar, one.shape), one, two, total)
    draws = np.stack(tuple(term[:, curves.rows].T for term in terms), axis=-1)
    best = np.stack(tuple(term[:, curves.best_row] for term in terms), axis=-1)
    return radius, draws, best


def plot_profiles(curves, observed_radii, stellar_mass):
    """Prepare all best-fit curves and the posterior bands for their sum."""
    radius, stellar, one, two, total = curve_terms(curves, observed_radii, stellar_mass)
    row = curves.best_row
    best = np.column_stack((stellar[:, 0], one[:, row], two[:, row], total[:, row]))
    quantiles = np.percentile(total[:, curves.rows], PERCENTILES, axis=1)
    return Profiles(radius, best, quantiles[:, :, None])


def split_one_halo(profiles, parameters, omega_m):
    """Split the best one-halo curve into its weighted centered and off-center terms."""
    log_mass, concentration, f_cen = np.atleast_2d(parameters).T[:3, :, None]
    mass = 10.0**log_mass
    r180 = halo_radius(mass, RHO_CRITICAL * omega_m)
    scale_radius, normalization = nfw_scales(mass, concentration, r180)
    centered = (normalization * (f_cen * nfw_esd(profiles.radius/scale_radius)))[0]
    stellar, one, two, total = profiles.best.T
    best = np.column_stack((stellar, centered, one - centered, two, total))
    return Profiles(profiles.radius, best, profiles.quantiles)


def read_profiles(path, fit_id=None):
    with Path(path).open() as stream:
        if fit_id is not None and stream.readline().strip() != f"# fit_id={fit_id}":
            raise ValueError("Fit summary and profile table have different fit IDs")
        table = np.loadtxt(stream, ndmin=2)
    expected = 1 + len(COMPONENTS)*len(COLUMNS)
    if table.shape[0] < 2 or table.shape[1] != expected:
        raise ValueError(f"Invalid ESD profile table: {path}")
    radius = table[:, 0]
    if not np.all(np.isfinite(table)) or radius[0] <= 0 or np.any(np.diff(radius) <= 0):
        raise ValueError(f"Invalid values or radius ordering in ESD profile table: {path}")
    stats = table[:, 1:].reshape(len(table), len(COMPONENTS), len(COLUMNS))
    return Profiles(radius, stats[:, :, 0], stats[:, :, 1:].transpose(2, 0, 1))
