# Configuration

Pass the configuration file path as the sole argument to `esd_bcg_mpi`. With no argument, the program looks for `esd.ini` in the launch directory and uses built-in defaults if the file is absent. Omitted keys use their defaults. Relative paths in the INI start from the launch directory.

## Run paths

| Key under `[run]` | What to set |
| --- | --- |
| `group-list`, `shear-list` | Paths to the two input lists. The defaults are `group_catalogs.list` and `shear_catalogs.list` in the launch directory. Both files are required. [Input formats](INPUTS.md) |
| `output` | Filename for the measured ESD table. The default is `ESD_results.dat`. Create its parent directory and use a distinct name for each result you want to keep. |
| `cache` | Optional path to save and reuse resolved centers. LWC and NWC add `_lwc` or `_nwc` before the extension. Create the cache's parent directory. |
| `scratch-dir` | Directory for temporary files. The default is `.bcg_scratch` beside the ESD table. On a multi-node run, use a shared location with sufficient free space. |
| `memory-mb` | Buffer budget per process in MiB. Default `512`, minimum `64`. |

## Group sample and measurement

| Key under `[measurement]` | Default | What to set |
| --- | --- | --- |
| `center` | `bcg` | `bcg` for the brightest member, `lwc` for the luminosity-weighted center, or `nwc` for the number-weighted center. |
| `z-min`, `z-max` | `0.3`, `0.4` | Group redshift range. The lower limit is included and the upper limit is excluded. |
| `logm-min`, `logm-max` | `13.0`, `14.2` | Group catalog mass range in base-10 log of solar masses divided by `h`. The lower limit is included and the upper limit is excluded. |
| `r-min`, `r-max` | `0.1`, `100.0` | Projected radius range in comoving Mpc divided by `h`. The lower limit is included and the upper limit is excluded. |
| `n-r-bins` | `12` | Number of radial bins from `6` to `256`. |
| `dz-min` | `0.2` | A source must be more than this redshift difference behind its group. |
| `n-jack` | `400` | Target number of sky regions for measurement uncertainties, at least `7`. A fit needs more usable regions than fitted radial bins. |
| `pdf-bins` | `4` | Keep the default for a standard run. |

Use matching cosmology settings when comparing samples. The defaults are `omega-m` at `0.315`, `h` at `0.674`, `sigma8` at `0.811`, `omega-b-h2` at `0.02237`, `spectral-index` at `0.965`, and `neutrino-mass` at `0.06` eV. The program generates its fit inputs from these settings.

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
| `walkers` | `72` | Even number, at least `10`. |
| `steps`, `burn-in` | `10000`, `1000` | Steps per walker and initial steps to discard. Use fewer discard steps than total steps. |
| `ensembles`, `seed` | `4`, `12771` | Number of independent fit runs and the nonnegative seed. |
| `profile-draws` | `4096` | Maximum number of fitted profile draws to save. The full parameter chain is saved separately. |
| `python-path` | unset | Optional directory containing Python packages for the fit. |

In `[plots]`, set `title` to override the default center-based title or set `esd-plot` to choose a different `.png` path. The default figure name follows the ESD table name. [Result files](RESULTS.md)
