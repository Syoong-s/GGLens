# Results

With `[run] output-dir = results` and `output = ESD_results.dat`, a successful `all` run produces these files.

| File | Contents |
| --- | --- |
| `measurement/ESD_results_measurement.bin` | All measured bins, both sets of jackknife replicas, lens-stack and calibration statistics, and cosmology |
| `fit/ESD_results/ESD_results_chain.bin` | All retained five-parameter samples and log probabilities in their original ensemble, step, and walker order, plus sampling diagnostics and the measurement reference |
| `fit/ESD_results/ESD_results_profile_draws.bin` | Joint ESD curve predictions at the selected chain positions and at the global best fit |
| `ESD_results_fit.png` | Tangential ESD used in the fit, best-fit stellar, centered one-halo, off-center one-halo, two-halo and total curves, and the total's 68% and 95% posterior bands |
| `ESD_results.dat` | All measured radial bins and measurement statistics |
| `ESD_results_posterior.png` | Five marginal posteriors with the 16th, 50th and 84th percentiles, and pairwise joint distributions |

`measure` writes the measurement snapshot. `fit` writes the two fit archives. `pull` reads those saved products and generates the three files in the result directory. The current fit directory is a relative symlink to the committed pair, so copy the fit directory with its symlink target when moving a run.

The binary snapshots preserve double precision. The measurement includes the source count after cuts, the number of all-radius sources used for self-calibration, the effective calibration coefficients, and the BCG photometry selection count when available. The total number of lens–source pairs is the sum of the per-bin counts. The jackknife replicas supply the fit covariance and both sets of ESD errors.

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
