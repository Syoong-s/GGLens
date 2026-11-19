# Input catalogs

Set `group-list` and `shear-list` in the `[run]` section. Use plain-text lists and whitespace-separated catalog records.

## DR9 Group, Galaxy, and Member catalogs

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

The program accept input in DR9 format (See [DR9 Format](dr9.md)). The tables below list the columns required by the program. Group IDs are local to a block. Galaxy ID `0` is valid.

| Catalog | Column, starting at 1 | Required value |
| --- | ---: | --- |
| Group | 1 | Group ID |
| Group | 5 | Group redshift |
| Group | 6 | Base-10 logarithm of the catalog halo mass in solar masses divided by `h` |
| Member | 1 | Galaxy ID from the Galaxy catalog |
| Member | 2 | Group ID from the Group catalog |
| Member | 3 | Brightest-member flag. BCG selection uses the member marked `1`. |

| Galaxy column, starting at 0 | Required value | Center |
| ---: | --- | --- |
| 0 | Galaxy ID | All |
| 1, 2 | Right ascension and declination in degrees | All |
| 5 | Positive, K-corrected `z`-band luminosity | LWC |
| 9 | Base-10 logarithm of stellar mass in solar masses | BCG |
| 14 | Photometry selection code | BCG |

Keep all member rows for LWC and NWC. For BCG, each selected group needs exactly one brightest-member row. Every required member must have a matching Galaxy ID.

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

Use degrees for RA and Dec, with RA in `[0, 360)` and Dec in `[-90, 90]`. The default photo-z cuts keep `0 < zp < 5`. For the header below, set `ra = 1`, `dec = 2`, `zp = 3`, `mag = 4`, and `expo-chi2 = -1` under `[external_columns]`. These are one-based positions, with `-1` selecting the last column for exposure quality.

```text
# ra dec zphot mag g1 g2 de h1 h2 gf1 gf2 half_light_flux half_light_area nstar imax jmax expo_chi2
```

Add the columns needed by any optional cuts you enable. Some use named columns, while magnitude, morphology, exposure ID, CCD ID, and exposure quality use column positions. [Source-cut columns](CONFIGURATION.md#source-cuts)
