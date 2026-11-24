# Results

With `[run] output-dir = results` and `output = ESD_results.dat`, a successful `all` run produces these files.

| File | Contents |
| --- | --- |
| `measurement/ESD_results_measurement.bin` | All measured bins, both sets of jackknife replicas, lens-stack and calibration statistics, and cosmology |
| `measurement/ESD_results_fd_points.dat` | FD test points for both Cartesian components when self-calibration is enabled |
| `fit/ESD_results/ESD_results_chain.bin` | All retained five-parameter samples and log probabilities in their original ensemble, step, and walker order, plus sampling diagnostics and the measurement reference |
| `fit/ESD_results/ESD_results_profile_draws.bin` | Joint ESD curve predictions at the selected chain positions and at the global best fit |
| `ESD_results_fit.png` | Tangential ESD used in the fit, best-fit stellar, centered one-halo, off-center one-halo, two-halo and total curves, and the total's 68% and 95% posterior bands |
| `ESD_results.dat` | All measured radial bins and measurement statistics |
| `ESD_results_posterior.png` | Five marginal posteriors with the 16th, 50th and 84th percentiles, and pairwise joint distributions |
| `ESD_results_fd.png` | FD test for both Cartesian components when the saved measurement used self-calibration and its point table is available |

`measure` writes the measurement snapshot. `fit` writes the two fit archives. `pull` reads the saved products and exports the plots and ESD table in the result directory. The current fit directory is a relative symlink to the committed pair, so copy the fit directory with its symlink target when moving a run.

The binary snapshots preserve double precision. The measurement includes the source count after cuts, the number of all-radius sources used for self-calibration, the effective calibration coefficients, and the BCG photometry selection count when available. The total number of lens–source pairs is the sum of the per-bin counts. The jackknife replicas supply the fit covariance and both sets of ESD errors.

## FD test points

With `[calibration] self-calibrate = true`, `measure` and `all` save the FD points as soon as the bin calculations finish. The table is available if the subsequent line fit or calibration validation fails.

| Columns | Contents |
| --- | --- |
| `component`, `bin` | Component 1 for `g1/gf1` or 2 for `g2/gf2`, and the bin number starting at 1 |
| `gf_low`, `gf_high`, `gf_mean` | Lower and upper FD bin edges and the mean FD of the sources in the bin |
| `shear`, `sigma` | PDF_SYM shear estimate and its curvature-based error bar |
| `n_sources` | Number of selected source rows in the bin, counting each row once per component across all matched lenses |
| `status` | `0` for a valid estimate, `1` for too few sources, or `2` for an invalid PDF_SYM fit |

An empty bin has `nan` for `gf_mean`. Bins with status `1` or `2` have `nan` for `shear` and `sigma`.

During `pull`, the FD plot uses the valid bins, their saved error bars, the calibration lines from the measurement snapshot, and a unit-slope reference. Keep the `*_fd_points.dat` table beside its matching `*_measurement.bin` when moving or reusing a measurement. A missing or mismatched table leaves the ESD and posterior exports available.

## ESD data table

Each row describes one radial interval. Radii are in comoving Mpc divided by `h`. ESD and its errors are in `h` solar masses per square parsec.

| Columns, starting at 1 | Contents |
| --- | --- |
| 1–3 | Lower edge, upper edge, and effective radius |
| 4 | Number of lens–source pairs |
| 5–6 | Tangential ESD and its error |
| 7–8 | Cross ESD and its error |

An unavailable value appears as `nan`. The table retains bins without complete jackknife measurements alongside the fitted bins.

## Reusing earlier results

A legacy `*_fit_input.dat` can still supply a fit. A later `pull` also needs its full companion ESD table to export every measured bin. Older results with an ESD table and corresponding `*_cov.dat` are supported with the original cosmology in the INI. Existing legacy fits with a matching `*_fit.json`, `*_profiles.dat`, and `*_chain.npz` can also be pulled.
