# Model, sampling, and products

## Model definition

The fit uses the BCG, LWC, and NWC conventions of Katsianis et al. The total predicted ESD is $\Delta\Sigma_\star+\Delta\Sigma_{\rm 1h}+\Delta\Sigma_{\rm 2h}$. The lens sample is fixed by the catalog redshift and mass cuts. The model mass is a fit parameter. [Katsianis et al., Section III and Appendix A](https://arxiv.org/abs/2207.12771)

The fitting code works in comoving `Mpc/h` and `h Msun/pc²`. The mean matter density in these coordinate and mass units is $\bar\rho_{m,0}=2.77536627\times10^{11}\Omega_m$. The point-mass term for BCG uses the linear mean stellar mass of participating BCGs. The stack converts the input masses from `Msun` to numeric `Msun/h` by multiplying by `h`. Evaluate $\Delta\Sigma_\star=M_\star/(10^{12}\pi R^2)$. LWC and NWC pass zero for $M_\star$.

### One halo

For $M=10^{\log_{10}M_{180m}}$, the code uses $r_{180m}=[3M/(4\pi\,180\bar\rho_{m,0})]^{1/3}$, $r_s=r_{180m}/c$, and $x=R/r_s$. The surface-density and excess-density kernels share the normalization $M/[2\pi r_s^2(\log(1+c)-c/(1+c))]$. The conversion to `h Msun/pc²` divides by $10^{12}$.

`nfw_projection` evaluates the analytic projected NFW surface-density kernel $F(x)$, with a series about $x=1$. `nfw_esd` computes the dimensionless excess as $2[\log(x/2)+H(x)]/x^2-F(x)$, using an algebraic rearrangement at small $x$ to limit cancellation.

The one-halo mixture has weight $f_{\rm cen}$ on the centered profile and $1-f_{\rm cen}$ on the offset profile. The offset distribution is Rayleigh with $R_{\rm sig}=\eta r_{180m}$. In units of $r_s$, the scale is $\sigma=\eta c$. Convolve the projected surface density with a two-dimensional isotropic Gaussian of that width. The angular integral uses a scaled modified Bessel function of order zero. Use 96-point Gauss–Legendre integration, then integrate the smoothed surface density inside the aperture before subtracting the local value. The offset ESD is tabulated at 128 logarithmic offsets and 640 logarithmic radii and interpolated with a bivariate spline. Use the centered limit below the tabulated offset range. [NFW kernels](../code/script/esd_nfw.py) · [Model assembly](../code/script/esd_fit.py)

### Two halos

At the mean group redshift of the participating lenses, rank 0 generates the nonlinear matter spectrum with CAMB's original Halofit. The cold matter density is $\Omega_mh^2-\Omega_bh^2-\sum m_\nu/93.14$ in the code's convention. Rescale the scalar amplitude by the squared ratio of the configured to the initial $\sigma_8(0)$ and recompute the spectrum from the transfer functions. Perform the spherical-Bessel transform with a 4,096-point, tapered logarithmic FFTLog grid. This gives $\xi_{mm}(r)$ in comoving coordinates.

The scale-dependent factor is $\zeta_0=(1+1.17\xi_{mm})^{1.49}/(1+0.69\xi_{mm})^{2.09}$. Freeze $\zeta$ below the last downward crossing of $\log_{10}(\zeta_0\xi_{mm})=0.9$. The broadcast array contains radii and $q(r)=\zeta(r)\xi_{mm}(r)$. The fit excludes separations inside $r_{180m}$ in the two-halo term.

For each projected radius $R$, set $u=\operatorname{asinh}(\Pi/R)$ and $u_{\min}=\operatorname{arcosh}(\max(1,r_{180m}/R))$. The projection combines an inner spherical integral with the line-of-sight tail.

$$
\Delta\Sigma_{\rm 2h}(R)=\frac{b_h\bar\rho_{m,0}}{10^{12}}
\left[
\frac{4}{R^2}\int_{r_{180m}}^R q(r)r^2\,dr
-2R\int_{u_{\min}}^{u_{\max}} q(R\cosh u)\cosh u\,e^{-2u}\,du
\right].
$$

The first integral contributes only for $R>r_{180m}$. The code uses PCHIP for the spherical primitive and a 513-point cumulative Simpson table per radius for the tail. Each proposal changes the mass-dependent lower limits and the linear amplitude $b_h$. [Correlation and projection](../code/script/esd_two_halo.py)

## Likelihood and chain

The likelihood is Gaussian with the complete covariance submatrix. `FitProblem` takes its Cholesky factor $L$ and scores a proposal using $-\lVert L^{-1}(x-m)\rVert^2/2$. The five sampling coordinates have the following open, flat prior ranges.

| Parameter | Range |
| --- | --- |
| $\log_{10}(M_{180m}/[h^{-1}M_\odot])$ | `(12, 16)` |
| $c$ | `(1, 20)` |
| $f_{\rm cen}$ | `(0, 1)` |
| $\eta=R_{\rm sig}/r_{180m}$ | `(0, 1)` |
| $b_h$ | `(0.3, 10)` |

`fit_esd` assigns contiguous ensemble indexes to MPI ranks. Each assigned ensemble initializes its NumPy generator from the configured seed and ensemble index and evaluates the log probability in batches through emcee. The bridge gathers the post-burn samples and log probabilities to rank 0. The summary reports acceptance and autocorrelation estimates per ensemble, including the retained-length-to-50-autocorrelation-times check. The best fit is the retained draw with maximum log probability, equivalent to minimum chi-squared under the flat priors. The reported degrees of freedom are the number of fitted radii minus five. [MPI fit](../code/src/esd_fit_mpi.cpp) · [Sampler](../code/script/esd_fit.py)

## Output data and rendering

`*_chain.npz` stores `posterior` with axes `(ensemble, retained step, walker, parameter)` and the corresponding `log_probability`, parameter `names`, observed `radii`, `measured`, `covariance`, zero-based `fit_bins`, and `center`. The JSON parameter percentiles use the full retained chain.

The profile grid is the union of the fitted effective radii and 180 logarithmically spaced plot points. Select up to `profile-draws` rows from the flattened chain without replacement, using the configured seed and a separate seed sequence. The `chain_indices` preserve the mapping to the full chain. Compute each component and their joint sum for every chosen draw. Profile percentiles are pointwise over the selected draws and use the sum for each draw when forming the total bands. The best curve has its own evaluation from the global best draw. The draw file uses `(draw, radius, component)` with component order `stellar`, `one_halo`, `two_halo`, `total`. The fit JSON also records best values at the observed radii.

The fit ID hashes the retained chain and the measured and plotted data that identify the fit. The profile table, profile archive, chain, and summary carry the same ID. The writers stage each file on its destination filesystem and publish the fit summary after the plot. The saved-figure reader validates the source IDs, observed radii, and best curve before rendering. It uses the same rendering routine and numeric values as the main run. [Curve export](../code/script/esd_profiles.py) · [Saved-figure reader](../code/script/esd_plot.py)

`multi_esd.py` compares historical 11-column and current 8-column ESD tables using the measured tangential values in column 5. It divides corresponding rows by the first table, takes error magnitudes from column 6, and reads radius from column 3. The reference is held fixed for this display. Use comparable radial rows for scientific comparisons. The vertical display is an `asinh` scale centered on ratio 1, and horizontal display offsets are independent of the input values. [Comparison script](../code/script/multi_esd.py)
