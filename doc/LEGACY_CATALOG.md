# Legacy i-band catalog layout

Source list: `/lustre/home/acct-phyzj/share/DES/catalogs/cat_iband_ori.list`.
Observed catalog: `/lustre/home/acct-phyzj/share/DES/catalogs/cat_iband_ori/subcat_001682.cat`.
Inspection was limited to its first header and first numeric row; both contain **44 fields**.

| 1-based column | Header name |
| ---: | --- |
| 1 | `flags_footprint` |
| 2 | `flags_foreground` |
| 3 | `flags_gold` |
| 4 | `ext_mash` |
| 5 | `ra` |
| 6 | `dec` |
| 7 | `bdf_mag_g` |
| 8 | `bdf_mag_err_g` |
| 9 | `bdf_mag_r` |
| 10 | `bdf_mag_err_r` |
| 11 | `bdf_mag_i` |
| 12 | `bdf_mag_err_i` |
| 13 | `bdf_mag_z` |
| 14 | `bdf_mag_err_z` |
| 15 | `bdf_mag_y` |
| 16 | `bdf_mag_err_y` |
| 17 | `dnf_z` |
| 18 | `dnf_zsigma` |
| 19 | `ccD_NUM` |
| 20 | `poly_chi2` |
| 21 | `xc` |
| 22 | `yc` |
| 23 | `sigma` |
| 24 | `nstar` |
| 25 | `imax` |
| 26 | `jmax` |
| 27 | `half_light_flux` |
| 28 | `half_light_area` |
| 29 | `flag` |
| 30 | `psf_FWHM` |
| 31 | `SNR_F` |
| 32 | `ra` |
| 33 | `dec` |
| 34 | `gf1` |
| 35 | `gf2` |
| 36 | `g1` |
| 37 | `g2` |
| 38 | `de` |
| 39 | `h1` |
| 40 | `h2` |
| 41 | `cos2` |
| 42 | `sin2` |
| 43 | `parity` |
| 44 | `Chi2` |

Rank 0 counts this catalog's 44 columns and uses the first layout for the whole run. The compiled positions select columns 5–6 for coordinates, 17 for photometric redshift, 11 for magnitude and the last column for exposure chi2. These positions can be overridden under `[external_columns]`. The magnitude and exposure chi2 columns are read when their cuts need them.

Global star-bar cuts use the selected magnitude, flux/area columns 27/28 and PSF size column 30. SNR is flux divided by the square root of area. Per-exposure star-bar fitting requires a configured exposure ID column. This catalog provides CCD identifiers.
