# Configuration

Pass the configuration file path to `esd_bcg_mpi`. When the path is omitted, the program looks for `esd.ini` in the launch directory and uses built-in defaults if the file is absent. Omitted keys use their defaults. Relative paths in the INI start from the launch directory, with the output path rules below.

## Run paths

| Key under `[run]` | What to set |
| --- | --- |
| `stage` | `all` runs `measure`, `fit`, and `pull` in order and is the default. `measure` creates the binary measurement snapshot. `fit` reads it and writes the full chain and profile draws. `pull` reads both stages and exports the figures and table. Use `--stage` to override the INI for one launch. |
| `group-list`, `shear-list` | Paths to the two input lists, required for `measure` and `all`. The defaults are `group_catalogs.list` and `shear_catalogs.list` in the launch directory. [Input formats](INPUTS.md) |
| `output-dir` | Optional result directory. The measurement snapshot goes in `measurement/`, the two fit archives go in `fit/<prefix>/`, and the three display products go in the result directory. |
| `output` | Name of the exported ESD table and source of the shared file prefix, by default `ESD_results.dat`. Use a filename with `output-dir`. With `output-dir` omitted, the parent of `output` serves as the result directory. |
| `cache` | Optional path to save and reuse resolved centers. LWC and NWC add `_lwc` or `_nwc` before the extension. |
| `scratch-dir` | Directory for temporary files. The default is `.bcg_scratch` inside `output-dir` or beside the ESD table when `output-dir` is omitted. On a multi-node run, use a shared location with sufficient free space. |
| `memory-mb` | Buffer budget per process in MiB. Default `512`, minimum `64`. |

## Group sample and measurement

| Key under `[measurement]` | Default | What to set |
| --- | --- | --- |
| `center` | `bcg` | `bcg` for the brightest member, `lwc` for the luminosity-weighted center, or `nwc` for the number-weighted center. |
| `z-min`, `z-max` | `0.3`, `0.4` | Group redshift range. The lower limit is included and the upper limit is excluded. |
| `logm-min`, `logm-max` | `13.0`, `14.2` | Group catalog mass range in base-10 log of solar masses divided by `h`. The lower limit is included and the upper limit is excluded. |
| `r-min`, `r-max` | `0.1`, `100.0` | Projected radius range in comoving Mpc divided by `h`. The lower limit is included and the upper limit is excluded. |
| `n-r-bins` | `12` | Number of radial bins, up to `256`. Measurement alone needs at least one bin. The fit needs at least six bins with complete jackknife measurements. |
| `dz-min` | `0.2` | A source must be more than this redshift difference behind its group. |
| `n-jack` | `400` | Target number of sky regions for measurement uncertainties. Use at least `2` for measurement only or `7` for `all`. The fit needs more usable regions than fitted radial bins. |
| `pdf-bins` | `4` | Keep the default for a standard run. |

Use matching cosmology settings when comparing samples. The defaults are `omega-m` at `0.315`, `h` at `0.674`, `sigma8` at `0.811`, `omega-b-h2` at `0.02237`, `spectral-index` at `0.965`, and `neutrino-mass` at `0.06` eV. The program generates its fit inputs from these settings.

## Group, Galaxy, and Member columns

Positions are one-based in all three sections. Negative positions count backward, with `-1` for the last column. The defaults match the DR9 layout. A position is checked at startup when the selected center uses it.

| Section | Key | Default | Value in the column |
| --- | --- | ---: | --- |
| `[group_columns]` | `id` | `1` | Positive Group ID. Set `0` to use the data row number. |
| `[group_columns]` | `redshift` | `5` | Group redshift. |
| `[group_columns]` | `log-halo-mass` | `6` | Group halo mass in the representation set by `mass-format`. |
| `[galaxy_columns]` | `id` | `1` | Nonnegative Galaxy ID. Set `0` to use the data row number. |
| `[galaxy_columns]` | `ra`, `dec` | `2`, `3` | Galaxy coordinates in degrees. |
| `[galaxy_columns]` | `luminosity` | `6` | Positive luminosity for LWC. |
| `[galaxy_columns]` | `log-stellar-mass` | `10` | BCG stellar mass in the representation set by `mass-format`. |
| `[galaxy_columns]` | `photometry-select` | `15` | BCG photometry selection code for the run log. Set `0` to omit the statistic. |
| `[member_columns]` | `galaxy-id`, `group-id` | `1`, `2` | IDs used to join each Member row to the Galaxy and Group tables. |
| `[member_columns]` | `brightest-flag` | `3` | BCG flag column. A selected group needs exactly one row whose flag equals `brightest-value`. |

Each section accepts `delimiter = whitespace` (default), `comma`, or `tab`. Set `header = true` to skip the first physical line of every file of that type. The default is `false`. Blank lines and lines beginning with `#` after optional whitespace are skipped in either case.

`[group_columns]` and `[galaxy_columns]` accept `mass-format = log10` (default) or `linear`. Group masses are in solar masses divided by `h` and enter the sample cuts as `log10` values. Galaxy stellar masses are in physical solar masses and enter the BCG stack as linear values. `[member_columns]` accepts an integer `brightest-value`, default `1`.

For `id = 0`, numbering starts at `1` in each file and counts every data record, including rows that fail a selection cut. The Member table uses these IDs in its configured join columns. [Input layouts and example](INPUTS.md#group-galaxy-and-member-catalogs)

## Shear column positions

Set each applicable key in `[external_columns]` to its one-based position in the shear catalogs.

| Key | Default | When to set it |
| --- | ---: | --- |
| `ra`, `dec`, `zp` | `5`, `6`, `17` | In every run, match your RA, Dec, and photometric redshift columns. |
| `mag` | `11` | Select the magnitude band for magnitude or star-bar cuts. |
| `ext-mash` | `4` | Set when using `ext-cut`. |
| `expo-num` | `0` (unset) | Supply the exposure ID position when using `expo-bar-mltp`. |
| `ccd-num` | `0` (unset) | Supply the `CCD_NUM` column position when using `ccd-cut`. |
| `expo-chi2` | `-1` | Supply the exposure quality position when using `chi2-thresh`. Use `-1` for the last column. |

## Source cuts

Set a threshold **below `-999.0`** to disable its cut. Values at or above `-999.0` are active, including `0` and negative thresholds. Omitted keys keep the defaults shown below. Each active cut needs its listed columns, including any columns shared with a disabled cut.

| Cut key | Default | Behavior | Required column or setting |
| --- | --- | --- | --- |
| `snrfcut` | `-9999.0` (off) | Keep values at least this threshold. | `snr_f` |
| `snrlow`, `snrhigh` | `10`, `-9999.0` (off) | Keep this source signal-to-noise range, including its endpoints. | `half_light_flux`, `half_light_area` |
| `zplow`, `zphigh` | `0`, `5` | Keep source photo-z strictly between the active limits. | `zp` |
| `mag-low`, `mag-high` | `10`, `30` | Keep the magnitude range, including its endpoints. | `mag` |
| `ext-cut` | `-9999.0` (off) | Keep the specified morphology value within `0.001`. | `ext-mash` |
| `starcut` | `16` | Require at least this number of stars. | `nstar` |
| `chi2-thresh`, `psf-chi2-mltp` | `0.01`, `-9999.0` (off) | Keep the corresponding quality value at or below the threshold. | `expo-chi2` and `poly_chi2`, respectively |
| `flagcut` | `-9999.0` (off) | Keep `flag × sigma` at or below the threshold. | `flag`, `sigma` |
| `imaxcut`, `jmaxcut` | `64`, `64` | Keep each value strictly below its active limit. | `imax`, `jmax` |
| `r-half-thresh` | `-9999.0` (off) | Require a source-to-PSF size ratio above this threshold. | `half_light_area`, `psf_fwhm` |
| `star-bar-mltp`, `expo-bar-mltp` | Both `-9999.0` (off) | Set the global or per-exposure star rejection cut. Per-exposure selection takes precedence when both are enabled. | `mag`, `half_light_flux`, `half_light_area`, `psf_fwhm`. The per-exposure cut also requires `expo-num`. |
| `ccd-cut` | `-9999.0` (off) | Exclude rows with any of the listed CCD IDs, across all exposures and shear catalogs. | `ccd-num` |
| `max-dup` | `5` | Set a positive integer maximum for adjacent rows with the same source position. Blocks over the limit are dropped. | `ra`, `dec` |

Enter any number of nonnegative CCD IDs on one line, separated by commas. For example, `ccd-cut = 2, 31, 53, 61` excludes those four CCDs from measurement and self-calibration. Set `ccd-num = 19` if your `CCD_NUM` is in the 19th column. A blank `ccd-cut` or a single value below `-999.0` disables it.

## Calibration

| Key under `[calibration]` | Default | What to set |
| --- | --- | --- |
| `fd-corrected` | `false` | Set to `true` if field distortion was already removed from input `g1` and `g2`. |
| `self-calibrate` | `true` | Estimate calibration values from this run. Set to `false` to use manual values. |
| `m1`, `m2`, `c1`, `c2` | `0` | Supply manual values when automatic calibration is off. Use `m1` and `m2` greater than `-1`. |
| `gf-lim` | `0.0015` | Set the maximum absolute field distortion for each component. A positive value also sets the automatic calibration range from `-gf-lim` to `+gf-lim`. With this cut disabled, automatic calibration uses `0.0015` for the half-range. |
| `selfcal-fd-bins` | `21` | Number of field-distortion intervals per component when automatic calibration is on. |
| `selfcal-pdf-bins`, `selfcal-min-count`, `selfcal-search-width` | `4`, `6`, `0.2` | Keep the defaults for a standard automatic calibration. |

## Fit and figure

| Key under `[fit]` | Default | What to set |
| --- | --- | --- |
| `input` | unset | In `stage = fit`, set a path to a `*_measurement.bin` from another run. The default is the snapshot derived from `[run] output` and `output-dir`. A later `pull` follows the fit's saved reference. Set this key in `pull` to locate the same snapshot after moving it. |
| `walkers` | `72` | Even number, at least `10`. |
| `steps`, `burn-in` | `10000`, `1000` | Steps per walker and initial steps to discard. Use fewer discard steps than total steps. |
| `ensembles`, `seed` | `4`, `12771` | Number of independent fit runs and the nonnegative seed. |
| `profile-draws` | `4096` | Maximum number of retained chain positions sampled to compute profile bands. The archive preserves every selected position and includes the global best curve. |
| `python-path` | unset | Optional directory containing Python packages for fitting or plotting. |

The saved measurement supplies the center, lens statistics, and cosmology. Earlier `*_fit_input.dat` snapshots and matching eight-column ESD tables with `*_cov.dat` can also be passed as `input`. For older ESD/covariance pairs, supply the original cosmology in `[measurement]`. To export every measured point after fitting from a legacy `*_fit_input.dat`, place its original ESD table beside it.

In `[plots]`, set `title` to override the center-based ESD figure title, or set `esd-plot` to choose a different `.png` path for `pull`. With `output-dir` set, a relative `esd-plot` starts from that directory. The default figure names follow the ESD table name. [Result files](RESULTS.md)
