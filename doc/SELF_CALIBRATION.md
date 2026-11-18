# Cartesian all-R FD self-calibration

Set `self-calibrate = true` under `[calibration]` in the INI to fit the Cartesian shear response to FD. Set `fd-corrected = true` when the input catalog already removed the FD signal. The two settings are independent.

## Restore raw moments

For each accepted source, define `B1 = de - h1` and `B2 = de + h1`. When `fd-corrected` is enabled, the code restores the moments before writing the source spool:

```text
G1_raw = G1_catalog + gf1 * B1 - gf2 * h2
G2_raw = G2_catalog + gf2 * B2 - gf1 * h2
```

This is the inverse of the catalog FD correction. The source keeps its original `de`, `h1`, `h2`, `gf1` and `gf2`. With raw input, the catalog `g1/g2` are used directly.

## Sample and calibration order

1. Stream the catalog and apply the quality and MAX_DUP cuts. Fit and apply the configured star-bar cuts.
2. For self-calibration, select the union of rows matching at least one lens in the configured radial and redshift ranges. A row matched to multiple lenses enters the union once. Distinct exposures remain separate measurement rows.
3. Measure raw Cartesian G1/G2 in uniform component FD bins, using responses `B1/B2`. The original occupancy refinement, threshold-20 search shrinking and 200-point quadratic PDF fit supply each bin's value and curvature error.
4. Fit each component using `g_measured = (1+m)*g_FD + c` and inverse-variance weights. Broadcast the four fitted coefficients. With self-calibration disabled, use the configured manual coefficients.
5. Apply the gains and intercepts to the original or restored moments, subtract the full FD response, and rotate each pair for ESD.

```text
G1_corr = (G1_raw - c1 * B1) / (1 + m1) - gf1 * B1 + gf2 * h2
G2_corr = (G2_raw - c2 * B2) / (1 + m2) - gf2 * B2 + gf1 * h2
```

The ESD PDF_SYM and jackknife errors condition on the effective fixed coefficients.

## Controls

| INI key under `[calibration]` | Default | Meaning |
| --- | ---: | --- |
| `fd-corrected` | false | Restore raw G1/G2 from FD-corrected catalog moments |
| `self-calibrate` | false | Fit both Cartesian components on the all-R union |
| `m1`, `m2`, `c1`, `c2` | 0 | Manual calibration for runs with self-calibration disabled |
| `gf-lim` | -9999.0 (source cut disabled) | Source FD cutoff and half-width of the FD grid. The example INI uses 0.0015. |
| `selfcal-fd-bins` | 21 | Count of equal-width bins for each of `gf1` and `gf2` |
| `selfcal-pdf-bins` | 4 | Number of magnitude cells in each FD-bin PDF |
| `selfcal-min-count` | 6 | Minimum global rows for a component FD bin |
| `selfcal-search-width` | 0.2 | Initial zero-centered shear search half-width |

For positive `gf-lim`, each component divides `[-gf-lim,+gf-lim)` into `selfcal-fd-bins` equal-width bins of width `2 * gf-lim / selfcal-fd-bins`. Each bin includes its left edge. If `gf-lim` is below `-999.0`, the source FD cut is disabled and the fit grid uses a half-width of `0.0015`. The fit needs at least two valid, distinct bins per component and positive gains. FD bin scratch is processed one cell at a time, within the configured per-rank buffer budget.

Rank 0 retains the FD measurements and fits in memory and renders the diagnostic with the ESD figure at the end. The FD figure defaults to `ESD_results_FD.png` beside `ESD_results.dat`. Set `[plots] fd-plot` to change this path. The ESD table records the effective coefficients and the number of selected FD rows.
