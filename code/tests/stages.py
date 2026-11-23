"""Compare a combined run with two ways to resume an existing ESD measurement."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


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


def check_fit(reference, other):
    expected_summary = json.loads(Path(f"{reference}_fit.json").read_text())
    assert json.loads(Path(f"{other}_fit.json").read_text()) == expected_summary
    assert expected_summary["cosmology"]["h"] == 0.7
    with np.load(f"{reference}_chain.npz") as expected, np.load(f"{other}_chain.npz") as actual:
        np.testing.assert_array_equal(actual["fit_bins"], [0, 1, 3, 4, 5, 6])
        for name in expected.files:
            np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)
    assert other.with_name(other.name + "_fit.png").is_file()


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

        write_inputs(root / "inputs")
        run("combined", "[run]\noutput-dir = combined\noutput = sample.dat\ncache = combined/cache/lenses.dat\n"
            + CATALOG_SETTINGS + FIT_SETTINGS)
        measured = run("measurement", "[run]\nstage = measure\noutput-dir = measure only\noutput = sample.dat\n"
                       + CATALOG_SETTINGS + "[fit]\nwalkers = 1\n[plots]\nesd-plot = unused.svg\n")
        assert "Measuring bin" in measured and "CAMB" not in measured
        directory = root / "measure only" / "measurement"
        products = {path.name: path.read_bytes() for path in directory.iterdir()}
        assert set(products) == {"sample.dat", "sample_cov.dat", "sample_jackknife.dat", "sample_fit_input.dat"}
        assert not (directory.parent / "sample_fit.json").exists()
        small_measurement = CATALOG_SETTINGS.replace("n-r-bins = 7\nn-jack = 12", "n-r-bins = 1\nn-jack = 2")
        run("small", "[run]\nstage = measure\noutput-dir = small\n" + small_measurement)
        assert (root / "small" / "measurement" / "ESD_results_fit_input.dat").is_file()
        (root / "inputs").rename(root / "offline inputs")
        (root / "saved measurement.dat").write_bytes(products["sample_fit_input.dat"])

        resumed = run("fit", "[run]\nstage = measure\noutput-dir = refit\noutput = sample.dat\n"
                      "group-list =\nshear-list =\n[measurement]\ncenter = nwc\nh = 0.5\nr-min = -1\n"
                      + FIT_SETTINGS + "input = saved measurement.dat\n", "--stage", "fit")
        assert "Reading ESD fit input" in resumed and "Measuring bin" not in resumed
        assert not (root / "refit" / "measurement").exists()
        reference = root / "combined" / "sample"
        check_fit(reference, root / "refit" / "sample")
        assert products == {path.name: path.read_bytes() for path in directory.iterdir()}

        legacy = root / "legacy"
        legacy.mkdir()
        legacy_table = [line for line in products["sample.dat"].decode().splitlines() if not line.startswith("# fit_input=")]
        (legacy / "sample.dat").write_text("\n".join(legacy_table) + "\n")
        (legacy / "sample_cov.dat").write_bytes(products["sample_cov.dat"])
        run("legacy", "[run]\nstage = fit\noutput = legacy/sample.dat\n[measurement]\n" + COSMOLOGY + FIT_SETTINGS)
        check_fit(reference, legacy / "sample")

        bad = root / "bad" / "measurement"
        bad.mkdir(parents=True)
        (bad / "ESD_results_fit_input.dat").write_bytes(products["sample_fit_input.dat"].replace(
            b"FQ_GGLENS_FIT_INPUT 1", b"FQ_GGLENS_FIT_INPUT 2", 1
        ))
        run("bad", "[run]\nstage = fit\noutput-dir = bad\n" + FIT_SETTINGS,
            error="Unsupported ESD fit input version")

    print("Combined, split and legacy ESD fit checks passed")


if __name__ == "__main__":
    main()
