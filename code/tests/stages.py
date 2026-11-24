"""Check independent stages, exact archive recovery and legacy fit inputs."""

import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))

from esd_archive import CurveArchive, load_fit, save_fit
from esd_plot import posterior_percentiles
from esd_profiles import plot_profiles, restore_draws, split_one_halo, stellar_profile


COSMOLOGY = """\
omega-m = 0.28
h = 0.7
"""

CATALOG_SETTINGS = f"""\
group-list = inputs/groups.list
shear-list = inputs/shear.list
memory-mb = 64
[measurement]
r-min = 0.2
r-max = 25.6
n-r-bins = 7
n-jack = 12
{COSMOLOGY}[external_columns]
ra = 1
dec = 2
zp = 3
[calibration]
self-calibrate = false
[cuts]
snrlow = -9999
starcut = -9999
chi2-thresh = -9999
imaxcut = -9999
jmaxcut = -9999
mag-low = -9999
mag-high = -9999
"""

FIT_SETTINGS = f"""\
[fit]
walkers = 12
steps = 12
burn-in = 4
ensembles = 2
profile-draws = 16
python-path = {os.environ.get('FQ_GGLENS_TEST_PYTHON_PATH', '')}
"""


def write_inputs(path):
    path.mkdir()
    (path / "groups.list").write_text("# group_cat galaxy_cat member_cat\n1\ngroups.dat\ngalaxies.dat\nmembers.dat\n")
    (path / "shear.list").write_text("sources0.dat\nsources1.dat\n")
    groups, galaxies, members = [], [], []
    header = "# ra dec zp g1 g2 de h1 h2 gf1 gf2"
    sources = [[header], [header]]
    rng = np.random.default_rng(2485)
    for identifier, ra in enumerate(np.arange(10, 360, 22.5), 1):
        groups.append(f"{identifier} 1 0 0 0.35 13.5")
        galaxies.append(f"{identifier} {ra} 0 0 0 1 0 0 0 11 0 0 0 0 1")
        members.append(f"{identifier} {identifier} 1")
        for radial_bin in (0, 1, 3, 4, 5, 6):
            radius = 0.2 * 2**(radial_bin + 0.5)
            offsets = rng.uniform(0.9, 1.1, 96)
            noise1, noise2 = rng.normal(0, 0.04, (2, 96))
            for offset, g1, g2 in zip(offsets, noise1, noise2):
                source_ra = (ra + np.rad2deg(offset * radius / 1000)) % 360
                sources[identifier % 2].append(f"{source_ra:.16g} 0 0.9 {g1:.16g} {g2:.16g} 1 0 0 0 0")
    for name, rows in zip(("groups", "galaxies", "members", "sources0", "sources1"),
                          (groups, galaxies, members, *sources)):
        (path / f"{name}.dat").write_text("\n".join(rows) + "\n")


def measurement_path(directory, name="sample"):
    return directory / "measurement" / f"{name}_measurement.bin"


def load_measurement(path):
    data = path.read_bytes()
    header = struct.Struct("<8s4I8d3q3I5d")
    (magic, version, center, n, regions, z, stellar_mass, *fields) = header.unpack_from(data)
    assert (magic, version, center) == (b"FQESDMEA", 1, 0)
    groups = np.frombuffer(data, dtype="<u4", count=regions, offset=header.size)
    offset = header.size + groups.nbytes
    edges = np.frombuffer(data, dtype="<f8", count=n + 1, offset=offset)
    offset += edges.nbytes
    row_type = np.dtype([("radius", "<f8"), ("count", "<i8"), ("tangential", "<f8"), ("cross", "<f8")])
    rows = np.frombuffer(data, dtype=row_type, count=n, offset=offset)
    replicas = np.frombuffer(data, dtype="<f8", count=2*n*regions, offset=offset + rows.nbytes)
    assert len(data) == offset + rows.nbytes + replicas.nbytes + 8
    return groups, z, stellar_mass, fields, edges, rows, replicas.reshape(2, n, regions)


def check_fd_points(root, run):
    inputs = root / "inputs"
    lens_ras = np.loadtxt(inputs / "galaxies.dat", usecols=1)
    radial_bins = (0, 1, 3, 4, 5, 6)
    counts = np.array([2400, 0, 2400, 4, 2400])
    edges = np.linspace(-0.0015, 0.0015, len(counts) + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    expected_counts = np.concatenate((counts, counts[::-1]))
    valid = expected_counts >= 6
    settings = CATALOG_SETTINGS.replace("inputs/shear.list", "inputs/fd.list").replace(
        "self-calibrate = false", f"self-calibrate = true\nselfcal-fd-bins = {len(counts)}"
    )
    for name, gain2 in (("fd-valid", 0.9), ("fd-rejected", -0.9)):
        (inputs / "fd.list").write_text(f"{name}.dat\n")
        rng = np.random.default_rng(3817)
        with (inputs / f"{name}.dat").open("w") as stream:
            stream.write("# ra dec zp g1 g2 de h1 h2 gf1 gf2\n")
            for bin_index, (field, count) in enumerate(zip(centers, counts)):
                noise = rng.normal(0, 0.025, (count // 2, 2))
                for row, (noise1, noise2) in enumerate(np.concatenate((noise, -noise))):
                    radial_bin = radial_bins[(row // len(lens_ras)) % len(radial_bins)]
                    radius = 0.2 * 2**(radial_bin + 0.5)
                    offset = 0.91 + 0.015 * bin_index + 1e-5 * row
                    ra = lens_ras[row % len(lens_ras)] + np.rad2deg(offset * radius / 1000)
                    g1 = 1.1 * field + 0.0001 + noise1
                    g2 = -gain2 * field - 0.0002 + noise2
                    stream.write(f"{ra:.17g} 0 0.9 {g1:.17g} {g2:.17g} 1 0 0 {field:.17g} {-field:.17g}\n")
        snapshot = measurement_path(root / name)
        error = "m1,m2 > -1" if gain2 < 0 else None
        run(name, f"[run]\nstage = measure\noutput-dir = {name}\noutput = sample.dat\n" + settings, error=error)
        points = np.genfromtxt(snapshot.with_name("sample_fd_points.dat"), names=True)
        np.testing.assert_array_equal(points["component"], np.repeat([1, 2], len(counts)))
        np.testing.assert_array_equal(points["bin"], np.tile(np.arange(1, len(counts) + 1), 2))
        np.testing.assert_array_equal(points["n_sources"], expected_counts)
        np.testing.assert_array_equal(points["status"], ~valid)
        np.testing.assert_allclose(points["gf_mean"], np.where(expected_counts, np.tile(centers, 2), np.nan),
                                   atol=1e-15)
        assert np.all(np.isfinite(points["shear"][valid])) and np.all(points["sigma"][valid] > 0)
        assert np.all(np.isnan(points["shear"][~valid])) and np.all(np.isnan(points["sigma"][~valid]))
        assert snapshot.is_file() == (error is None)


def check_fd_pull(root, run):
    source = measurement_path(root / "fd-valid")
    points = source.with_name("sample_fd_points.dat")
    original_points = points.read_bytes()
    settings = ("[run]\noutput-dir = fd refit\noutput = refitted.dat\n"
                "[calibration]\nself-calibrate = false\n" + FIT_SETTINGS)
    run("fd-fit", settings + f"input = {source}\n", "--stage", "fit")
    output = root / "fd refit" / "refitted_fd.png"
    assert not output.exists()
    run("fd-pull", settings, "--stage", "pull")
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert points.read_bytes() == original_points
    output.unlink()

    rejected = root / "fd-rejected" / "measurement" / points.name
    points.write_bytes(rejected.read_bytes())
    run("fd-pull-mismatched", settings, "--stage", "pull")
    assert not output.exists()
    points.unlink()
    log = run("fd-pull-missing", settings, "--stage", "pull")
    assert "missing point table" in log
    assert {p.name for p in output.parent.iterdir() if p.is_file()} == {
        "refitted.dat", "refitted_fit.png", "refitted_posterior.png"
    }
    points.write_bytes(original_points)


def check_fit(reference, other, same_measurement=True):
    expected, expected_curves = load_fit(reference / "fit" / "sample")
    actual, curves = load_fit(other / "fit" / "sample")
    if same_measurement:
        assert actual.fit_id == expected.fit_id
    assert actual.shape == expected.shape and actual.seed == expected.seed
    np.testing.assert_array_equal(actual.expand(), expected.expand())
    np.testing.assert_array_equal(actual.diagnostics, expected.diagnostics)
    np.testing.assert_array_equal(curves.extra_radii, expected_curves.extra_radii)
    np.testing.assert_array_equal(curves.chain_indices, expected_curves.chain_indices)
    np.testing.assert_array_equal(curves.rows, expected_curves.rows)
    np.testing.assert_array_equal(curves.components, expected_curves.components)
    assert curves.best_row == expected_curves.best_row


def running_sum(values):
    value = 0.0
    for item in values:
        value += item
    return value


def legacy_inputs(directory, table, measurement):
    groups, z, mass, fields, _, rows, replicas = measurement
    regions = len(groups)
    cosmology = fields[:6]
    bins = [i for i, row in enumerate(rows) if np.isfinite(row["tangential"]) and np.all(np.isfinite(replicas[0, i]))]
    covariance = np.full((len(rows), len(rows)), np.nan)
    centered = {i: replicas[0, i] - running_sum(replicas[0, i]) / regions for i in bins}
    for i in bins:
        for j in bins:
            covariance[i, j] = running_sum(centered[i] * centered[j]) * (regions - 1.0) / regions
    directory.mkdir()
    header = f"center=bcg n_jack={regions} fit_bins_zero_based=" + " ".join(map(str, bins))
    np.savetxt(directory / "sample_cov.dat", covariance, fmt="%.17g", header=header)
    legacy = "\n".join(line for line in table.splitlines() if not line.startswith("# measurement=")) + "\n"
    (directory / "sample.dat").write_text(legacy)
    subset = covariance[np.ix_(bins, bins)]
    with (directory / "sample_fit_input.dat").open("w") as stream:
        stream.write(f"FQ_GGLENS_FIT_INPUT 1\ncenter bcg\nsample {sum(groups)} {regions} {z:.17g} {mass:.17g}\n")
        stream.write("cosmology " + " ".join(f"{value:.17g}" for value in cosmology) + "\n")
        stream.write(f"radial_bins {len(rows)}\nfit_bins {len(bins)}\nindex_radius_esd\n")
        for i in bins:
            stream.write(f'{i} {rows[i]["radius"]:.17g} {rows[i]["tangential"]:.17g}\n')
        stream.write("covariance\n")
        np.savetxt(stream, subset, fmt="%.17g")
    return bins


def check_lossless_codec(root):
    chain = np.arange(144, dtype=np.float64).reshape((1, 4, 6, 6))
    chain[:, 1] = chain[:, 0]
    chain[:, 3, ::2] = chain[:, 2, ::2]
    chain[0, 1, 0, 0] = -0.0
    chain[0, 2, 0, -1] = 999.0
    chain[0, 1, 1, -1] = 999.0
    diagnostics = np.array([[.4, 1., 2., 3., 4., np.nan]])
    draws = CurveArchive(2, np.empty(0), np.array([1, 2, 3], dtype="<u4"),
                         np.array([0, 0, 1], dtype="<u4"), 2, np.arange(12., dtype="<f8").reshape(2, 2, 3))
    source = root / "measurement.bin"
    archive = root / "fit" / "codec"
    context = {"measurement_id": 123, "measurement_path": source, "title": "", "burn_in": 12, "seed": 42}
    save_fit(archive, chain, diagnostics, draws, context)
    saved, curves = load_fit(archive)
    np.testing.assert_array_equal(saved.expand().view(np.uint64), chain.view(np.uint64))
    np.testing.assert_array_equal(saved.diagnostics, diagnostics)
    assert saved.best() == (7, 999.0)
    index, parameters = saved.best_sample()
    assert index == 7
    np.testing.assert_array_equal(parameters, chain.reshape(-1, 6)[index])
    for i in range(5):
        np.testing.assert_array_equal(posterior_percentiles(saved.states[i], saved.weights),
                                      np.percentile(chain[..., i], [16, 50, 84]))
    radius, complete, best = restore_draws(curves, [1., 2.], 4.)
    np.testing.assert_array_equal(complete[0], complete[1])
    np.testing.assert_array_equal(complete[:, :, 0], np.broadcast_to(stellar_profile(radius, 4.), (3, 2)))
    np.testing.assert_array_equal(complete[:, :, 3], (complete[:, :, 0] + complete[:, :, 1]) + complete[:, :, 2])
    profiles = plot_profiles(curves, radius, 4.)
    np.testing.assert_array_equal(profiles.best, best)
    np.testing.assert_array_equal(profiles.quantiles[:, :, 0],
                                  np.percentile(complete[:, :, 3], [2.5, 16, 50, 84, 97.5], axis=0))
    parameters = [13.6, 6.5, 0.5, 0.18, 2.0]
    split = split_one_halo(profiles, parameters, 0.28)
    np.testing.assert_array_equal(split.best[:, [0, 3, 4]], best[:, [0, 2, 3]])
    np.testing.assert_allclose(split.best[:, 1] + split.best[:, 2], best[:, 1])
    np.testing.assert_array_equal(split.quantiles, profiles.quantiles)
    parameters[2] = 0.0
    split = split_one_halo(profiles, parameters, 0.28)
    np.testing.assert_array_equal(split.best[:, 1], np.zeros_like(radius))
    np.testing.assert_array_equal(split.best[:, 2], best[:, 1])
    previous = archive.resolve()
    save_fit(archive, chain, diagnostics, draws, context)
    assert not previous.exists()
    np.testing.assert_array_equal(saved.expand().view(np.uint64), chain.view(np.uint64))
    curve_path = archive / "codec_profile_draws.bin"
    blob = curve_path.read_bytes()
    curve_path.write_bytes(blob[:-1] + bytes([blob[-1] ^ 1]))
    try:
        load_fit(archive)
    except ValueError as error:
        assert "checksum mismatch" in str(error)
    else:
        raise AssertionError("Corrupted curve archive was accepted")


def main():
    executable = Path(sys.argv[1]).resolve()
    command = [*sys.argv[2:], "-n", "2", executable]
    with tempfile.TemporaryDirectory(prefix="fq-esd-stages-") as tmp:
        root = Path(tmp)
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                   MPLCONFIGDIR=str(root / "matplotlib"))

        def run(name, settings, *options, error=None):
            ini = root / f"{name}.ini"
            ini.write_text(settings)
            result = subprocess.run([*command, *options, ini], cwd=root, env=env,
                                    capture_output=True, text=True, timeout=120)
            if error is not None:
                assert result.returncode != 0 and error in result.stderr, result.stdout + result.stderr
            else:
                assert result.returncode == 0, result.stdout + result.stderr
            return result.stdout

        check_lossless_codec(root)
        write_inputs(root / "inputs")
        check_fd_points(root, run)
        combined = root / "combined"
        run("combined", "[run]\noutput-dir = combined\noutput = sample.dat\ncache = combined/cache/lenses.dat\n"
            + CATALOG_SETTINGS + FIT_SETTINGS)
        displays = {"sample.dat", "sample_fit.png", "sample_posterior.png"}
        assert {path.name for path in combined.iterdir() if path.is_file()} == displays
        assert {path.name for path in (combined / "fit" / "sample").iterdir()} == {
            "sample_chain.bin", "sample_profile_draws.bin"
        }

        measured = run("measurement", "[run]\nstage = measure\noutput-dir = measure only\noutput = sample.dat\n"
                       + CATALOG_SETTINGS + "[fit]\nwalkers = 1\n[plots]\nesd-plot = unused.svg\n")
        assert "Measuring bin" in measured and "CAMB" not in measured
        measure_dir = root / "measure only"
        snapshot = measurement_path(measure_dir)
        assert {p.name for p in snapshot.parent.iterdir()} == {"sample_measurement.bin"}
        assert not any(p.is_file() for p in measure_dir.iterdir())
        saved_snapshot = snapshot.read_bytes()
        assert saved_snapshot == measurement_path(combined).read_bytes()
        measurement = load_measurement(snapshot)
        groups, _, _, fields, edges, rows, replicas = measurement
        assert fields[1] == 0.7 and fields[6:9] == [16*6*96, -1, 16]
        assert sum(groups) == 16 and len(groups) == 12 and np.isnan(replicas[:, 2]).all()
        table = np.loadtxt(combined / "sample.dat")
        np.testing.assert_array_equal(table[:, 0], edges[:-1])
        np.testing.assert_array_equal(table[:, 1], edges[1:])
        np.testing.assert_array_equal(table[:, 2], rows["radius"])
        np.testing.assert_array_equal(table[:, 3], rows["count"])
        np.testing.assert_array_equal(table[:, 4], rows["tangential"])
        np.testing.assert_array_equal(table[:, 6], rows["cross"])

        small_measurement = CATALOG_SETTINGS.replace("n-r-bins = 7\nn-jack = 12", "n-r-bins = 1\nn-jack = 2")
        run("small", "[run]\nstage = measure\noutput = small/ESD_results.dat\n" + small_measurement)
        small_snapshot = measurement_path(root / "small", "ESD_results")
        assert small_snapshot.is_file()
        (root / "inputs").rename(root / "offline inputs")
        check_fd_pull(root, run)

        refit = root / "refit"
        refit_run = "[run]\nstage = measure\noutput-dir = refit\noutput = sample.dat\n"
        unrelated = "group-list =\nshear-list =\n[measurement]\ncenter = nwc\nh = 0.5\nr-min = -1\n"
        resumed = run("fit", refit_run + unrelated + FIT_SETTINGS + f"input = {snapshot}\n", "--stage", "fit")
        assert "Reading ESD fit input" in resumed and "Measuring bin" not in resumed
        assert not (refit / "measurement").exists() and not any(p.is_file() for p in refit.iterdir())
        check_fit(combined, refit)
        assert snapshot.read_bytes() == saved_snapshot

        pull_settings = refit_run + unrelated + FIT_SETTINGS.replace("walkers = 12", "walkers = 1")
        pulled = run("pull", pull_settings, "--stage", "pull")
        assert "Sources after cuts: 9216" in pulled and "select>=1: 16" in pulled
        assert "Measuring bin" not in pulled and "CAMB" not in pulled
        assert {path.name for path in refit.iterdir() if path.is_file()} == displays
        for name in displays:
            assert (refit / name).read_bytes() == (combined / name).read_bytes()
        run("pull-mixed", pull_settings + f"input = {small_snapshot}\n", "--stage", "pull",
            error="different measurement snapshot")

        legacy = root / "legacy"
        bins = legacy_inputs(legacy, (combined / "sample.dat").read_text(), measurement)
        assert bins == [0, 1, 3, 4, 5, 6]
        run("legacy-file", "[run]\nstage = fit\noutput = legacy/sample.dat\n[measurement]\n"
            + COSMOLOGY + FIT_SETTINGS + "input = legacy/sample.dat\n")
        check_fit(combined, legacy, same_measurement=False)
        stored, _ = load_fit(legacy / "fit" / "sample")
        np.testing.assert_array_equal(stored.legacy_cosmology, fields[:6])
        run("legacy-file-pull", "[run]\nstage = pull\noutput = legacy/sample.dat\n" + FIT_SETTINGS)
        assert (legacy / "sample_fit.png").read_bytes() == (combined / "sample_fit.png").read_bytes()
        run("legacy-snapshot", "[run]\nstage = fit\noutput = legacy/sample.dat\n" + FIT_SETTINGS)
        check_fit(combined, legacy, same_measurement=False)
        run("legacy-pull", "[run]\nstage = pull\noutput = legacy/sample.dat\n" + FIT_SETTINGS)
        np.testing.assert_array_equal(np.loadtxt(legacy / "sample.dat"), table)

        bad = root / "bad"
        bad_snapshot = measurement_path(bad, "ESD_results")
        bad_snapshot.parent.mkdir(parents=True)
        bad_snapshot.write_bytes(saved_snapshot[:-9] + bytes([saved_snapshot[-9] ^ 1]) + saved_snapshot[-8:])
        run("bad", "[run]\nstage = fit\noutput-dir = bad\n" + FIT_SETTINGS, error="archive checksum")

        bad_legacy = root / "bad legacy.dat"
        bad_legacy.write_text((legacy / "sample_fit_input.dat").read_text().replace(
            "FQ_GGLENS_FIT_INPUT 1", "FQ_GGLENS_FIT_INPUT 2", 1))
        run("bad-legacy", "[run]\nstage = fit\noutput-dir = bad\n" + FIT_SETTINGS + f"input = {bad_legacy}\n",
            error="Unsupported ESD fit input version")

    print("FD points and plots, combined and split stages, lossless archive, and legacy ESD checks passed")


if __name__ == "__main__":
    main()
