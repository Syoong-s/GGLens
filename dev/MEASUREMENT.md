# ESD and covariance

## Pair geometry and units

`Cosmology` tabulates the flat matter-plus-Lambda comoving distance in `Mpc/h`. It uses 4,001 points through redshift 5, Simpson integration per interval, and cubic Hermite interpolation with exact end slopes. The integrand is $2998/E(z)$, where $E(z)^2=\Omega_m(1+z)^3+1-\Omega_m$.

For each source chunk, build a KD tree of source sky vectors. Query a chord distance of $2\sin(R_{\max}/(2\chi_l))$ around each lens. Retain pairs with $z_s>z_l+\mathrm{dz\_min}$ and a projected radius in $[R_{\min},R_{\max})$. The radius uses the great-circle angle.

$$
R=\chi_l\,\operatorname{atan2}\left(\left\lVert\hat n_s\times\hat n_l\right\rVert,\;\hat n_s\cdot\hat n_l\right).
$$

For the source-local angle $\phi$ from east toward north, let $c_n=\cos(n\phi)$ and $s_n=\sin(n\phi)$. The stored moments and responses use these conventions.

$$
\begin{aligned}
G_t &= -(G_1c_2+G_2s_2), & G_\times &= G_1s_2-G_2c_2,\\
U_t &= -(h_1c_4+h_2s_4), & \beta &= \frac{10^{12}\chi_l(\chi_s-\chi_l)(1+z_l)}{1.661\times10^{18}\chi_s},\\
D_t &= \beta(de+U_t), & D_\times &= \beta(de-U_t).
\end{aligned}
$$

Here $\beta$ is inverse comoving critical density in `pc²/(h Msun)`. The source catalog photo-z supplies $z_s$. Distances and responses use the group redshift for all three center choices. The bin coordinates are logarithmically spaced and the effective radius is $\exp(\sum_i\log R_i/N_{\rm pair})$. [Distance and search code](../include/bcg_physics.hpp) · [Pair creation](../src/esd_pairs.cpp)

## Scratch and sky regions

The run creates one unique scratch subtree and an owned subdirectory per rank. Quality-selected source rows and optional star-bar metadata use aligned binary streams. After calibration, each radial bin has a pair stream. A record stores `gt`, `gx`, `det`, and `dex` as doubles, a 32-bit lens index, and padding, for 40 bytes total. Each bin tracks pair counts, summed log radius, and participating lens flags. The measurement loads a bin into memory when its records fit the resident quota, and streams it otherwise. Per-bin data and the owned run scratch are removed on successful completion.

Build the region list from the union of lenses in radial bins with enough global pairs. NGC and SGC contribute to one sky sample. Rank 0 seeds spherical k-means with a fixed-seed sample, up to the requested count or the number of eligible lenses. Worker ranks assign by maximum dot product and reduce the vector sums. Lloyd iterations stop when the labels stabilize or reach 32 iterations. Empty regions are removed and labels are compacted. Every radius uses the same group region labels. The fixed fit inputs use the same lens sample and average over lenses with equal weight. [Scratch code](../src/bcg_stream.cpp) · [Pair store](../src/bcg_pair_store.cpp) · [Region assignment](../src/esd_regions.cpp)

## PDF_SYM estimates

For each radial bin and component, use the residual $y=G-\Delta\Sigma D$. Select trial values of $\Delta\Sigma$ by minimizing the sign-count statistic

$$
\chi^2(\Delta\Sigma)=\sum_b
\frac{(N_{b,+}-N_{b,-})^2}{2(N_{b,+}+N_{b,-})}.
$$

Only populated bins contribute. Exact zero residuals are omitted in ESD trials. The initial magnitude edges are estimated separately for the full and deleted-region samples from a shared 2,048-cell log-$|G|$ histogram plus a zero cell. The full sample uses the sum over regions. A deletion subtracts the region's histogram and moments. Each of the configured `pdf-bins` edges targets a quantile $i/(\mathrm{pdf\_bins}+1)$ and the combined edges define `pdf-bins + 1` magnitude bins including overflow.

Trial sign counts are exact at the chosen edges. For each pair, the signed bin is monotone across ordered trials. The implementation counts intervals of trial indices using a difference array. It accumulates full-sample counts at the union of all edges and omitted counts at the corresponding region's edges. Prefix sums recover counts per trial and threshold, after which each replicate uses the full counts minus its omitted counts.

Each sample's initial center is $\sum DG/\sum D^2$. The bracket uses a scatter estimate from $\sum D^2$, $\sum DG$, and $\sum G^2$. The shared 41-point coarse grid can widen up to five times. A 61-point grid refines the range around the coarse minima. A local quadratic supplies an interior minimum when valid, with the grid minimum as fallback. Estimates with insufficient pairs, zero response, or a final minimum at the grid boundary remain unavailable. [Edge estimates](../src/bcg_jackknife_edges.cpp) · [Trial histograms](../src/bcg_jackknife_histogram.cpp) · [Trial search](../src/bcg_jackknife_fit.cpp)

## Full radial covariance

For $K$ regions, estimate the covariance from the full vector $x^{(-k)}$ of replicate tangential ESD values.

$$
C_{ij}=\frac{K-1}{K}\sum_{k=1}^{K}
\left(x_i^{(-k)}-\bar{x}_i\right)
\left(x_j^{(-k)}-\bar{x}_j\right).
$$

Every fitted bin requires a finite full value and all $K$ finite deletion values. Rows and columns unavailable for the fit are marked `nan` in the full matrix. The ESD table uses $\sqrt{C_{ii}}$ for the tangential error and a separate deletion variance for the cross error. The fit maps the valid radial indices to a dense covariance submatrix. It requires more regions than fitted bins, at least six bins, and a positive-definite submatrix. [Covariance and text output](../src/esd_covariance.cpp) · [Fit inputs](../src/esd_fit_mpi.cpp)
