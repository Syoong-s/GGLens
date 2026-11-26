"""Exercise the MPI lens loader with the same sample in several catalog layouts."""

import csv
import io
import math
from pathlib import Path
import subprocess
import sys
import tempfile


GROUPS = (
    ((1, 0.2, 13.5), (2, 0.35, 13), (3, 0.37, 15)),
    ((1, 0.4, 13.5), (2, 0.375, 14), (3, 0.31, 12)),
)

# ID, RA, Dec, luminosity, log stellar mass, photometry code.
GALAXIES = (
    (0, 359, 0, 1, 10, 0),
    (2**53 + 1, 1, 0, 3, 11, 1),
    (17, 42, 0, 5, 10, 2),
    (2**53, 98, 0, 2, 9, 0),
    (37, 102, 0, 4, 12, 3),
)
MEMBERS = ((1, 2), (4, 5))  # One-based galaxy rows with the BCG listed last.

WIDE = """[group_columns]
id = 3
redshift = -1
log-halo-mass = 1
[galaxy_columns]
id = -2
ra = -1
dec = 1
luminosity = 5
log-stellar-mass = 3
photometry-select = -3
delimiter = comma
header = true
[member_columns]
galaxy-id = -1
group-id = 3
brightest-flag = 1
brightest-value = 7
delimiter = comma
header = true
"""

IMPLICIT_GROUP = """[group_columns]
id = 0
redshift = 1
log-halo-mass = 2
mass-format = linear
delimiter = comma
header = true
"""

IMPLICIT = IMPLICIT_GROUP + """[galaxy_columns]
id = 0
ra = -1
dec = 1
luminosity = 4
log-stellar-mass = 3
photometry-select = 0
mass-format = linear
delimiter = tab
header = true
[member_columns]
brightest-value = 7
"""

MINIMAL = IMPLICIT_GROUP + """[galaxy_columns]
id = 0
ra = 1
dec = 2
"""


def write_table(path, rows, separator=" ", header=False):
    preamble = ["   # comment preceding the first data row", ""]
    if header:
        preamble.insert(0, separator.join(f"header_{column}" for column in range(len(rows[0]))))
    # Comments, nonuniform line lengths, and no final newline exercise byte partition alignment.
    between_rows = "\n\t\n# " + "padding " * 23 + "\n"
    lines = []
    for row in rows:
        if separator == ",":
            out = io.StringIO()
            csv.writer(out, quoting=csv.QUOTE_ALL, lineterminator="").writerow(row)
            lines.append(out.getvalue())
        else:
            lines.append(separator.join(map(str, row)))
    path.write_text("\n".join(preamble) + "\n" + between_rows.join(lines))


def write_sample(root, layout):
    root.mkdir()
    implicit = layout in {"implicit", "minimal"}
    wide = layout == "wide"
    galaxy_rows = []
    for identifier, ra, dec, lum, log_mass, select in GALAXIES:
        if wide:
            row = (dec, 'spare, "quoted"', log_mass, ra + 5, lum) + ("",) * 16 + (select, identifier, ra)
        elif layout == "implicit":
            row = (dec, "", 10**log_mass, lum, ra)
        elif layout == "minimal":
            row = (ra, dec)
        else:
            row = (identifier, ra, dec, 0, 0, lum, 0, 0, 0, log_mass, 0, 0, 0, 0, select)
        galaxy_rows.append(row)
    galaxy_separator = "," if wide else "\t" if layout == "implicit" else " "
    write_table(root / "galaxies", galaxy_rows, galaxy_separator, wide or layout == "implicit")

    manifest = ["# group_cat galaxy_cat member_cat"]
    for block, (groups, members) in enumerate(zip(GROUPS, MEMBERS), 1):
        group_path = f"groups{block}"
        member_path = f"members{block}"
        manifest.extend((str(block), group_path, "galaxies", member_path))
        group_rows = []
        for identifier, z, mass in groups:
            if wide:
                group_rows.append((mass, "unused", identifier, z))
            elif implicit:
                group_rows.append((z, 10**mass))
            else:
                group_rows.append((identifier, len(members), 0, 0, z, mass, 0))
        write_table(root / group_path, group_rows, "," if implicit else " ", implicit)

        member_rows = []
        for index, row_number in enumerate(members):
            identifier = row_number if implicit else GALAXIES[row_number - 1][0]
            if wide:
                member_rows.append((7 if index else 2, "", 2, identifier))
            elif layout == "minimal":
                member_rows.append((identifier, 2))
            else:
                central_flag = 7 if implicit else 1
                member_rows.append((identifier, 2, central_flag if index else 2, 0))
        write_table(root / member_path, member_rows, "," if wide else " ", wide)
    (root / "catalogs.list").write_text("\n".join(manifest) + "\n")


def load(root, launcher, ranks, center, settings, cache=False, error=None):
    cache_line = "cache = 'lenses #1.dat' # Keep the hash in the filename\n" if cache else ""
    (root / "run.ini").write_text(
        f'[run]\ngroup-list = "catalogs.list" # Lens inputs\n{cache_line}'
        f"[measurement]\ncenter = {center} # Group center\nn-r-bins = 12 # Radial bins\n{settings}"
    )
    result = subprocess.run(launcher(ranks), cwd=root, capture_output=True, text=True, timeout=30)
    if error is not None:
        assert result.returncode != 0 and error in result.stderr, result.stdout + result.stderr
        return
    assert result.returncode == 0, result.stdout + result.stderr
    lenses = []
    for line in result.stdout.splitlines():
        if line.startswith("lens "):
            _, block, group, ra, dec, z, mass, galaxy, stellar, select = line.split()
            lenses.append((int(block), int(group), float(ra), float(dec), float(z),
                           float(mass), int(galaxy), float(stellar), int(select)))
    return lenses, result.stdout


def check_centers(lenses, center, implicit=False, ra_offset=0):
    assert len(lenses) == len(GROUPS), lenses
    for block, (lens, groups, members) in enumerate(zip(lenses, GROUPS, MEMBERS), 1):
        _, z, mass = groups[1]
        ra, galaxy_id, stellar, select = 0, -1, 0, 0
        if center == "bcg":
            galaxy_id, ra, _, _, log_stellar, select = GALAXIES[members[-1] - 1]
            ra += ra_offset
            stellar = 10**log_stellar
            if implicit:
                galaxy_id, select = members[-1], 0
        else:
            angles = [math.radians(GALAXIES[row - 1][1]) for row in members]
            weights = [GALAXIES[row - 1][3] if center == "lwc" else 1 for row in members]
            x = sum(w * math.cos(a) for w, a in zip(weights, angles))
            y = sum(w * math.sin(a) for w, a in zip(weights, angles))
            ra = math.degrees(math.atan2(y, x)) % 360
        assert lens[:2] == (block, 2), lens
        assert abs(math.remainder(lens[2] - ra, 360)) < 1e-11, lens
        assert lens[3:] == (0, z, mass, galaxy_id, stellar, select), lens


def main():
    executable = Path(sys.argv[1]).resolve()

    def launcher(ranks):
        return [*sys.argv[2:], "-n", str(ranks), executable, "run.ini"]

    with tempfile.TemporaryDirectory(prefix="fq-catalogs-") as tmp:
        for layout, settings in (("legacy", ""), ("wide", WIDE), ("implicit", IMPLICIT)):
            root = Path(tmp) / layout
            write_sample(root, layout)
            for center in ("bcg", "lwc", "nwc"):
                for ranks in (1, 4):
                    lenses, _ = load(root, launcher, ranks, center, settings)
                    check_centers(lenses, center, implicit=layout == "implicit")

        root = Path(tmp) / "minimal"
        write_sample(root, "minimal")
        lenses, _ = load(root, launcher, 4, "nwc", MINIMAL)
        check_centers(lenses, "nwc", implicit=True)
        load(root, launcher, 4, "lwc", MINIMAL, error="galaxy_columns.luminosity")

        root = Path(tmp) / "wide"
        _, log = load(root, launcher, 1, "bcg", WIDE, cache=True)
        assert "from cache" not in log, log
        assert (root / "lenses #1.dat").is_file()
        lenses, log = load(root, launcher, 4, "bcg", WIDE, cache=True)
        assert "from cache" in log, log
        check_centers(lenses, "bcg")
        lenses, log = load(root, launcher, 4, "bcg", WIDE.replace("ra = -1", "ra = 4"), cache=True)
        assert "from cache" not in log, log
        check_centers(lenses, "bcg", ra_offset=5)
        load(root, launcher, 4, "bcg", WIDE.replace("ra = -1", "ra = -25"), error="galaxy_columns.ra")
    print("Catalog checks passed for DR9, reordered CSV, row IDs, minimal layouts, and cache changes.")


if __name__ == "__main__":
    main()
