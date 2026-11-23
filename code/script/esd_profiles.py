"""Posterior component curves and the table consumed by the fit plot."""

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

COMPONENTS = ("stellar", "one_halo", "two_halo", "total")
PERCENTILES = (2.5, 16, 50, 84, 97.5)
COLUMNS = ("best", "p2p5", "p16", "p50", "p84", "p97p5")


# ==========================================
# Class: Profiles
# Method: Keep one best curve and five marginal percentiles per radius and component.
# ==========================================
@dataclass(frozen=True)
class Profiles:
    radius: np.ndarray
    best: np.ndarray
    quantiles: np.ndarray


# ==========================================
# Function: write_atomic
# Method: Close a complete product in the output directory before publishing its final name.
# ==========================================
def write_atomic(path, mode, writer):
    path = Path(path)
    with TemporaryDirectory(dir=path.parent, prefix=f".{path.name}.") as directory:
        temp = Path(directory) / "product"
        with temp.open(mode) as stream:
            writer(stream)
        temp.replace(path)


# ==========================================
# Function: curve_array
# Method: Compute the three terms and their sum for every supplied joint parameter draw.
# ==========================================
def curve_array(model, theta, radius, two_halo):
    terms = model.components(theta, radius, two_halo)
    return np.stack((*terms, sum(terms)), axis=-1)


# ==========================================
# Function: evaluate_profiles
# Method: Evaluate a reproducible subset on a smooth grid that also includes all fitted radii.
# ==========================================
def evaluate_profiles(problem, posterior, best_index, limit, seed):
    from esd_two_halo import TwoHalo

    count = min(limit, len(posterior))
    if count == len(posterior):
        indices = np.arange(count, dtype=np.int64)
    else:
        rng = np.random.default_rng(np.random.SeedSequence([seed, 4796]))
        indices = np.sort(rng.choice(len(posterior), count, replace=False))
    smooth = np.geomspace(problem.r.min(), problem.r.max(), 180)
    radius = np.union1d(smooth, problem.r)
    model = problem.model
    two_halo = TwoHalo(radius, model.correlation, model.mean_density)
    best = curve_array(model, posterior[best_index], radius, two_halo)[0]

    draws = np.empty((count, len(radius), len(COMPONENTS)))
    batch = 256
    for start in range(0, count, batch):
        stop = min(count, start + batch)
        draws[start:stop] = curve_array(model, posterior[indices[start:stop]], radius, two_halo)
    if not (np.all(np.isfinite(best)) and np.all(np.isfinite(draws))):
        raise ValueError("Non-finite ESD profile in the retained posterior")
    profiles = Profiles(radius, best, np.percentile(draws, PERCENTILES, axis=0))
    return profiles, indices, draws


# ==========================================
# Function: write_profiles
# Method: Save the plotted curves and the corresponding joint draws with their chain indices.
# ==========================================
def write_profiles(table_path, draws_path, profiles, indices, draws, center, seed, fit_id):
    stats = np.concatenate((profiles.best[:, :, None], profiles.quantiles.transpose(1, 2, 0)), axis=2)
    columns = ["R_Mpc_over_h"] + [f"{term}_{stat}" for term in COMPONENTS for stat in COLUMNS]
    table = np.column_stack((profiles.radius, stats.reshape(len(profiles.radius), -1)))
    header = f"fit_id={fit_id}\n" + " ".join(columns) + f"\ncenter={center}  DeltaSigma in h M_sun/pc^2"
    write_atomic(table_path, "w", lambda stream: np.savetxt(stream, table, fmt="%.17g", header=header))
    write_atomic(draws_path, "wb", lambda stream: np.savez_compressed(
        stream, radius=profiles.radius, component_names=COMPONENTS, center=center,
        chain_indices=indices, draws=draws, seed=seed, fit_id=fit_id,
    ))


# ==========================================
# Function: read_profiles
# Method: Restore the published plot curves directly from the portable text table.
# ==========================================
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
