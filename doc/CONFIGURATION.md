# Configuration

Set values in `esd.ini`. Use `FQ_GGLENS_INI` to select another file. Omitted keys use their defaults. Relative paths in the INI start from the launch directory.

## Run paths

| Key under `[run]` | What to set |
| --- | --- |
| `group-list`, `shear-list` | Paths to the two [input lists](INPUTS.md). Both are required. |
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
| `r-min`, `r-max` | `0.1`, `5.0` | Projected radius range in comoving Mpc divided by `h`. The lower limit is included and the upper limit is excluded. |
| `n-r-bins` | `10` | Number of radial bins from `6` to `256`. |
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
| `expo-num` | unset | Supply the exposure ID position when using `expo-bar-mltp`. |
| `expo-chi2` | `-1` | Supply the exposure quality position when using `chi2-thresh`. Use `-1` for the last column. |

## Source cuts

Use the `[cuts]` section for optional source filters. Omit a cut or give it a value **below `-999.0`** to disable it. Values at or above `-999.0` are active, including `0` and negative thresholds. The program defaults leave these cuts disabled. The full example INI contains active example thresholds.

| Cut key | Behavior | Required column or setting |
| --- | --- | --- |
| `snrfcut` | Keep values at least this threshold. | `snr_f` |
| `snrlow`, `snrhigh` | Keep this source signal-to-noise range, including its endpoints. | `half_light_flux`, `half_light_area` |
| `zplow`, `zphigh` | Keep source photo-z strictly between the active limits. | `zp` |
| `mag-low`, `mag-high` | Keep the magnitude range, including its endpoints. | `mag` |
| `ext-cut` | Keep the specified morphology value within `0.001`. | `ext-mash` |
| `starcut` | Require at least this number of stars. | `nstar` |
| `chi2-thresh`, `psf-chi2-mltp` | Keep the corresponding quality value at or below the threshold. | `expo-chi2` and `poly_chi2`, respectively |
| `flagcut` | Keep `flag × sigma` at or below the threshold. | `flag`, `sigma` |
| `imaxcut`, `jmaxcut` | Keep each value strictly below its active limit. | `imax`, `jmax` |
| `r-half-thresh` | Require a source-to-PSF size ratio above this threshold. | `half_light_area`, `psf_fwhm` |
| `star-bar-mltp`, `expo-bar-mltp` | Set the global or per-exposure star rejection cut. Per-exposure selection takes precedence when both are enabled. | `mag`, `half_light_flux`, `half_light_area`, `psf_fwhm`. The per-exposure cut also requires `expo-num`. |
| `max-dup` | Set a positive integer maximum for adjacent rows with the same source position. Blocks over the limit are dropped. | `ra`, `dec` |

## Calibration

| Key under `[calibration]` | Default | What to set |
| --- | --- | --- |
| `fd-corrected` | `false` | Set to `true` if field distortion was already removed from input `g1` and `g2`. |
| `self-calibrate` | `false` | Set to `true` to estimate calibration values in this run. |
| `m1`, `m2`, `c1`, `c2` | `0` | Supply manual values when automatic calibration is off. Use `m1` and `m2` greater than `-1`. |
| `gf-lim` | disabled | Set the maximum absolute field distortion for each component. A positive value also sets the automatic calibration range from `-gf-lim` to `+gf-lim`. With this cut disabled, automatic calibration uses `0.0015` for the half-range. |
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
