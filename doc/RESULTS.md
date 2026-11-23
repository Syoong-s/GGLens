# Results

With `[run] output-dir = results` and `output = ESD_results.dat`, the program writes the following files under `results/`. Changing `output` changes the shared prefix. The figure path can also be set separately in `[plots]`.

| File | Contents |
| --- | --- |
| `measurement/ESD_results.dat` | Measured tangential and cross ESD, uncertainties, and pair counts by radius |
| `measurement/ESD_results_cov.dat` | Tangential ESD covariance in the same radial order as the measurement table. The header lists the zero-based rows used for the fit. |
| `measurement/ESD_results_jackknife.dat` | Tangential and cross ESD with each sky region omitted in turn |
| `measurement/ESD_results_fit_input.dat` | Self-contained input for a later five-parameter fit. It stores the fitted radial bins, tangential ESD, covariance submatrix, center, sample size, mean redshift, mean stellar mass, and cosmology. |
| `ESD_results_fit.json` | Center, sample size, cosmology, best-fit parameters, parameter intervals, and fit diagnostics |
| `ESD_results_chain.npz` | Retained parameter draws, measured values used in the fit, and their covariance |
| `ESD_results_profiles.dat` | Best-fit values and uncertainty bands for every fitted component and the total curve |
| `ESD_results_profile_draws.npz` | Saved draws of the individual curves. Draw indices identify the corresponding rows in the full parameter chain. |
| `ESD_results_fit.png` | Measured tangential ESD, total best-fit curve, and 68% and 95% bands |

`stage = measure` writes the files under `measurement/`. `stage = fit` reads the saved input and writes the fit products. Without `output-dir`, all files follow the parent directory of `output`.

## Measurement table

Each row describes one radial interval. Radii are in comoving Mpc divided by `h`. ESD and its errors are in `h` solar masses per square parsec.

| Columns, starting at 1 | Contents |
| --- | --- |
| 1–3 | Lower edge, upper edge, and effective radius |
| 4 | Number of lens–source pairs |
| 5–6 | Tangential ESD and its error |
| 7–8 | Cross ESD and its error |

An unavailable value appears as `nan`. Use the fit-bin indices in the covariance header to select the rows included in the fit.

## Plot-ready curves

The first column of `*_profiles.dat` is the plotted radius. Read each curve's six columns in the order `best`, `p2p5`, `p16`, `p50`, `p84`, `p97p5`.

| Column suffix | Meaning |
| --- | --- |
| `best` | Curve evaluated at the best-fit parameters |
| `p2p5`, `p97p5` | Lower and upper limits of the 95% band |
| `p16`, `p84` | Lower and upper limits of the 68% band |
| `p50` | Median of the saved profile draws |

The four curves are labeled `stellar`, `one_halo`, `two_halo`, and `total`. The stellar values are zero for LWC and NWC. Use the `total` columns for the total curve and bands. Keep the files with the same `fit_id` together when plotting or comparing a saved run.
