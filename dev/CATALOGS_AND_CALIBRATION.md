# Catalogs and calibration

## DR9 joins

The manifest scopes selected groups by `(block ID, Group ID)`. Rank 0 resolves catalog paths against the manifest directory, canonicalizes them, and sends the resolved manifest to the other ranks. The Group scan applies half-open redshift and catalog-mass limits. It retains Group ID, group redshift, and the selection mass. The selected groups are gathered before resolving coordinates.

`Dr9Rows` assigns a row to the rank whose byte interval contains the row's start. It keeps the first row after a newline and skips a leading fragment when the boundary cuts through a row. A rank can finish a row beyond its upper boundary. The reader skips blank lines and comments and handles an unterminated final row. IDs are parsed as integers from the records.

For BCG, each block matches exactly one Member row with brightest flag `1` per selected group. Requests are merged by canonical Galaxy path, including paths shared by NGC and SGC. Ranks read the corresponding explicit Galaxy IDs and gather the matching positions, photometric codes, and linear stellar masses computed as `10**log_stellar_mass`.

For LWC and NWC, every selected Member row contributes a Galaxy ID and the index of its group in the common lens table. For each unique Galaxy path, an exact shared ID selection uses a dense bitmap for suitable ID spans or a sparse set for wide spans. Requests are distributed by `Galaxy ID % number_of_ranks`. During the scan, workers route only matching positions to these owners, using batches of up to 262,144 data rows per reader. Owners detect missing or repeated IDs and accumulate group sums. MPI reduces the sums and every rank resolves the same centers.

For a member at right ascension $a$ and declination $d$, the sky vector is $\hat{n}=(\cos d\cos a,\cos d\sin a,\sin d)$. LWC uses the Galaxy catalog's column 5 luminosity as the positive weight in $\sum_i w_i\hat{n}_i$. NWC uses unit weights. Convert the summed direction to RA with `atan2(y,x)` and to Dec with `atan2(z,hypot(x,y))`. Wrap RA to `[0, 360)`. A group with no well-defined direction raises an error.

The `DR9_LENSES_V2` cache key includes the center, selection limits, block IDs, and input paths, sizes, and modification times. Rank 0 checks row counts and mode-specific data, then broadcasts a valid lens table. Cache misses write a complete replacement through a unique temporary directory on the same filesystem. [Manifest and cache](../src/bcg_dr9_manifest.cpp) · [Lens loading](../src/esd_lenses.cpp) · [Centroids](../src/bcg_dr9_centroids.cpp)

## Source schema and cuts

Rank 0 uses the first shear file's first line for the column count and mapping. It reserves the configured one-based external positions first, then looks backwards for each required name among unassigned columns. The first file sets the column width and order expected for the whole run. Readers skip their headers and apply the broadcast mapping. Non-finite selected values or those with absolute value above `1e300` cause the row to be dropped. Malformed numeric values in mapped columns raise an error.

Each source file is assigned by list position modulo the number of source workers. Rank 0 joins the source read for a single-rank run. The duplicate cut compares consecutive RA and Dec against the start of their current group to a tolerance of `1e-7` degree. Groups over the limit are dropped before per-row cuts. Duplicate grouping ends at each file boundary.

All optional cuts use `threshold >= -999.0` as the enable condition. Map each new named field in the first-header reader and tie its read dependency to the same cut predicate used by the row selector. The source quantities defined by the selector are

$$
\mathrm{SNR}=\frac{\mathrm{half\_light\_flux}}{\sqrt{A}},\qquad
d_{1/2}=2(0.2628)\sqrt{A/\pi},\qquad
s_{\rm rel}=\frac{d_{1/2}-\mathrm{PSF\ FWHM}}{\mathrm{PSF\ FWHM}},
$$

where $A$ is `half_light_area`. The relative size, selected magnitude, and SNR enter the star-bar stage. The exposure variant also uses a near-integer exposure ID between 1 and 10,000,000.

The older 44-column source layout placed external RA, Dec, photo-z, and i-band magnitude at columns 5, 6, 17, and 11. Its last column carried exposure chi-squared. The header rule is independent of width, and an enabled exposure star-bar cut needs the actual exposure ID position in the INI. [Header mapping](../src/bcg_shear_catalog.cpp) · [Cut logic](../src/bcg_shear_cuts.cpp)

## Star-bar selection

Fit the star bar after duplicate and simple row cuts. The size histogram uses 100 bins of width 0.04 over `[-2, 2)`. The magnitude histogram uses 20 unit bins starting at 10, with the reference's truncation convention. In the global mode, choose the strongest three-bin concentration in the permitted size range among magnitude bins with at least 100 rows. Collect nearby peaks within 0.05 in size and concentration of at least 0.6. Estimate the center and sample scatter using rows within 0.1 of the initial peak. The cut is the mean plus the configured multiplier times the scatter, with the code's fallback threshold when the peak or scatter estimate is insufficient. Global rejection removes rows at or below the cutoff when SNR is positive.

Per-exposure fitting uses its own size window around `[-0.5, 0.2]` when selecting peaks. The peak and mean stages use SNR above 40. Three clipping passes update the means and scatter, with a width floor of 0.005 on the scatter. The scatter stage uses positive SNR. Exposure IDs are collected with a bitmap, and the fits are evaluated in bounded batches. Failed exposures use the mean and scatter averaged over successful exposure fits, or the fixed fallback of 0.5 and 0.05 if every fit fails. The per-exposure cut takes precedence when enabled. It removes rows at or below the threshold for their exposure. [Star-bar code](../src/bcg_star_bar.cpp)

## Field distortion

Define $B_1=de-h_1$ and $B_2=de+h_1$. If `fd-corrected` is enabled, restore the raw moments before the source spool is written.

$$
G_{1,\rm raw}=G_{1,\rm cat}+gf_1B_1-gf_2h_2,\qquad
G_{2,\rm raw}=G_{2,\rm cat}+gf_2B_2-gf_1h_2.
$$

With self-calibration enabled, collect the union of post-cut source rows that make at least one pair in the selected radial and redshift ranges. Multiple lens matches of one row contribute once to the union. Distinct exposures remain separate rows. Each Cartesian component uses equal-width bins over `[-gf-lim, +gf-lim)`, with `selfcal-fd-bins` intervals. When the source FD cut is disabled, the grid half-width is 0.0015. Representative FD values are component-wise means within the retained calibration bins.

Measure each valid bin with the corresponding raw Cartesian moment and response $B_1$ or $B_2$. The FD solver uses exact local order statistics with bounded sorting and merged scratch runs. It averages eligible ranks' initial edges or gathers the sparse sample when every worker has too few rows. The global occupancy adjustment targets equal counts with a 0.03 fractional tolerance. The sign cost counts zero residuals in occupancy without adding to the sign difference. The search shrinks from a zero-centered interval using a cost difference of 20, then fits 200 points with a quadratic. A valid interior minimum with curvature $a>0$ has error $1/\sqrt{2a}$.

Rank 0 fits $g_{\rm measured}=(1+m)gf+c$ with inverse-variance weights and broadcasts the coefficients. Each component needs at least two valid, distinct FD points, and each gain $1+m$ must be positive. With self-calibration off, the INI supplies the coefficients. Apply them before pair rotation.

$$
\begin{aligned}
G_{1,\rm final} &= \frac{G_{1,\rm raw}-c_1B_1}{1+m_1}-gf_1B_1+gf_2h_2,\\
G_{2,\rm final} &= \frac{G_{2,\rm raw}-c_2B_2}{1+m_2}-gf_2B_2+gf_1h_2.
\end{aligned}
$$

The reported ESD errors condition on these fixed coefficients. The ESD header gives the effective coefficients and, when enabled, the calibration range and union size. [Restoration and correction](../include/bcg_calibration.hpp) · [FD line fit](../src/selfcalibration.cpp) · [FD sign solver](../src/selfcal_pdf.cpp)
