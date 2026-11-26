# Input catalogs

Set `group-list` and `shear-list` in the `[run]` section. Use plain-text lists and one record per line in the catalog files.

## Group, Galaxy, and Member catalogs

The group list starts with the header shown below. Each block has four lines in order. Give it a unique positive number, followed by the Group, Galaxy, and Member paths. Resolve relative paths from the list's directory. Multiple blocks can use the same Galaxy file.

```text
# group_cat galaxy_cat member_cat
1
catalogs/DESIDR9_NGC_group
catalogs/DESIDR9_galaxy
catalogs/iDESIDR9_NGC_1
2
catalogs/DESIDR9_SGC_group
catalogs/DESIDR9_galaxy
catalogs/iDESIDR9_SGC_1
```

Set positions in `[group_columns]`, `[galaxy_columns]`, and `[member_columns]`. Each catalog type can have its own column count, order, and delimiter. Rank 0 reads the first data row of the first file of each type to establish the three widths. Every file and data row of a given type must use that layout. [Column settings and DR9 defaults](CONFIGURATION.md#group-galaxy-and-member-columns) · [Original DR9 layout](dr9.md)

With `delimiter = comma`, a field can be double-quoted and can contain commas and escaped double quotes on a single line. Comma and tab separators preserve empty fields, so unused empty columns keep their positions. A required numeric field must contain a number.

Group IDs are local to the block, and explicit Galaxy ID `0` is valid. With `id = 0` in a Group or Galaxy section, the first data record of each file has ID `1`. The configured header, comments, and blank lines do not receive IDs. Member records use the resulting numbers to reference those files.

Keep all member rows for LWC and NWC. For BCG, each selected group needs exactly one brightest-member row. Every required member must have a matching Galaxy ID.

For a BCG run with headered CSV files, the following settings map a Group file with `redshift,halo_mass` and a Galaxy file with `ra,dec,stellar_mass`. Both store linear masses. The Member file has `galaxy_id,group_id,central`, with `7` marking its brightest members.

```ini
[group_columns]
id = 0
redshift = 1
log-halo-mass = 2
mass-format = linear
delimiter = comma
header = true

[galaxy_columns]
id = 0
ra = 1
dec = 2
log-stellar-mass = 3
photometry-select = 0
mass-format = linear
delimiter = comma
header = true

[member_columns]
galaxy-id = 1
group-id = 2
brightest-flag = 3
brightest-value = 7
delimiter = comma
header = true
```

## Shear catalogs

The shear list contains each single shear catalog path per line. Relative paths are resolved from the list's directory. Blank lines and lines beginning with `#` can be used in the list.

```text
shear/tile_001.cat
shear/tile_002.cat
```

Each catalog starts with a single whitespace-separated header, optionally beginning with `#`. Use the same column order for every shear file in a run. A standard single shear catalog contains one row per source with the following columns.

```text
# (External catalog fields) EXPO_NUM ccD_NUM poly_chi2 xc yc sigma nstar imax jmax half_light_flux half_light_area flag psf_FWHM SNR_F ra dec gf1 gf2 g1 g2 de h1 h2 cos2 sin2 parity Chi2
```

The required named shear columns for every run are **`g1`, `g2`, `de`, `h1`, `h2`, `gf1`, and `gf2`**. The default cuts also require `half_light_flux`, `half_light_area`, `nstar`, `imax`, and `jmax`, plus the configured magnitude and exposure quality columns.

Names are case-insensitive. If an imported column repeats a shear name, put the intended shear column later in the header.

Use degrees for RA and Dec, with RA in `[0, 360)` and Dec in `[-90, 90]`. The default photo-z cuts keep `0 < zp < 5`. For the header below, set `ra = 1`, `dec = 2`, `zp = 3`, and `expo-chi2 = -1` under `[external_columns]`. These are one-based positions, with `-1` selecting the last column for exposure quality. Set `mag-stellar-locus = 4` to fit the star bar with this magnitude. Use `ext-mag-i = [4, 10, 30]` for a magnitude cut on that column.

```text
# ra dec zphot mag g1 g2 de h1 h2 gf1 gf2 half_light_flux half_light_area nstar imax jmax expo_chi2
```

Add the columns needed by any optional cuts you enable. Magnitudes, magnitude errors, photo-z error, morphology, exposure ID, CCD ID, and exposure quality use column positions. Use `ext-` rules with one-based positions to filter on additional numeric columns. [Source-cut columns](CONFIGURATION.md#source-cuts) · [Additional column cuts](CONFIGURATION.md#additional-column-cuts)
