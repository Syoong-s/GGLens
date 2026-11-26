#pragma once
#include "bcg_calibration.hpp"
#include "bcg_dr9.hpp"
#include "bcg_shear_cuts.hpp"
#include <initializer_list>

namespace bcg_esd::defaults {
// INI discovery and paths
inline constexpr char ini_name[] = "esd.ini";
inline constexpr char group_list[] = "group_catalogs.list";
inline constexpr char shear_list[] = "shear_catalogs.list";
inline constexpr char output[] = "ESD_results.dat";
inline constexpr char cache[] = "";
inline constexpr char scratch_dir[] = "";
inline constexpr int memory_mib = 512;

// Lens sample, radial bins, cosmology, and PDF_SYM
inline constexpr bcg::Center center = bcg::Center::bcg;
inline constexpr double zmin = 0.3, zmax = 0.4, mmin = 13.0, mmax = 14.2;
inline constexpr double rmin = 0.1, rmax = 100.0, dzmin = 0.2, dz_err_k = 0.0;
inline constexpr double omega_m = 0.315, h = 0.674, sigma8 = 0.811;
inline constexpr double omega_b_h2 = 0.02237, spectral_index = 0.965, neutrino_mass = 0.06;
inline constexpr int n_r_bins = 12, n_jack = 400, pdf_bins = 4;

// DR9 positions are one-based for all three tables.
inline constexpr bcg::LensColumns lens_columns{
    {1, 5, 6, {}},
    {1, 2, 3, 6, 10, 15, {}},
    {1, 2, 3, {}}
};

// FD restoration and Cartesian calibration
inline constexpr bool fd_corrected = false, self_calibrate = true;
inline constexpr int selfcal_fd_bins = 21, selfcal_pdf_bins = 4, selfcal_min_count = 6;
inline constexpr double selfcal_search_width = 0.2;
inline constexpr double gf_lim = 0.0015;
inline constexpr double selfcal_fd_limit = 0.0015; // When the source FD cut is disabled.
inline constexpr bcg::ShearCalibration calibration{0.0, 0.0, 0.0, 0.0}; // m1, m2, c1, c2

// Independent emcee ensembles, with steps counted per walker.
inline constexpr int fit_walkers = 72, fit_steps = 10000, fit_burn = 1000;
inline constexpr int fit_ensembles = 4, fit_seed = 12771;
inline constexpr int fit_profile_draws = 4096;
inline constexpr char fit_python_path[] = "";
inline constexpr char fit_plot[] = "";
inline constexpr char plot_title[] = ""; // Derived from the center when empty.

// Positions are one-based. CatalogCombiner appends the exposure chi2.
inline constexpr auto external_columns = [] {
    bcg::ExternalColumns cols{};
    cols.ra = 5;
    cols.dec = 6;
    cols.zp = 17;
    cols.mag_stellar_locus = 11;
    cols.expo_num = 0; // Set explicitly for the per-exposure bar.
    cols.ccd_num = 0; // Set explicitly for the CCD cut.
    cols.expo_chi2 = -1; // Last column.
    return cols;
}();

// ==========================================
// Function: Default source cuts initializer
// Method: Set source selection defaults to the values in esd.example.ini.
// ==========================================
inline constexpr auto cuts = [] {
    bcg::CutThresholds c{};
    c.snrfcut = -9999.0;
    c.snrlow = 10.0;
    c.snrhigh = -9999.0;
    c.starcut = 16.0;
    c.chi2_thresh = 0.01;
    c.flagcut = -9999.0;
    c.zplow = 0.0;
    c.zphigh = 5.0;
    c.gf_lim = gf_lim;
    c.imaxcut = 64.0;
    c.jmaxcut = 64.0;
    c.psf_chi2_mltp = -9999.0;
    c.r_half_thresh = -9999.0;
    c.star_bar_mltp = -9999.0;
    c.expo_bar_mltp = -9999.0;
    c.max_dup = 5.0;
    return c;
}();

inline constexpr std::initializer_list<int> excluded_ccds = {}; // Equivalent to ccd-cut = [].
} // namespace bcg_esd::defaults
