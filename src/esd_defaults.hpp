#pragma once
#include "bcg_calibration.hpp"

namespace bcg_esd::defaults {
// INI discovery and paths
inline constexpr char ini_name[] = "esd.ini";
inline constexpr char ini_environment[] = "FQ_GGLENS_INI";
inline constexpr char group_list[] = "";
inline constexpr char shear_list[] = "";
inline constexpr char output[] = "ESD_results.dat";
inline constexpr char cache[] = "";
inline constexpr char scratch_dir[] = "";
inline constexpr int memory_mib = 512;

// Lens sample, radial bins, cosmology, and PDF_SYM
inline constexpr double zmin = 0.3, zmax = 0.4, mmin = 13.0, mmax = 14.2;
inline constexpr double rmin = 0.1, rmax = 5.0, dzmin = 0.2;
inline constexpr double omega_m = 0.315, h = 0.674, nfw_c = 7.44;
inline constexpr double conc_min = 1.0, conc_max = 40.0;
inline constexpr char conc_model[] = "duffy08_200m";
inline constexpr int n_r_bins = 10, n_jack = 100, pdf_bins = 4;

// FD restoration and Cartesian calibration
inline constexpr bool fd_corrected = false, self_calibrate = false;
inline constexpr int selfcal_fd_bins = 21, selfcal_pdf_bins = 4, selfcal_min_count = 6;
inline constexpr double selfcal_search_width = 0.2;
inline constexpr bcg::ShearCalibration calibration{0.0, 0.0, 0.0, 0.0}; // m1, m2, c1, c2

// Plot paths default to the result stem plus .png and _FD.png.
inline constexpr char esd_plot[] = "", fd_plot[] = "";
inline constexpr char plot_title[] = "BCG-centered group lensing";
inline constexpr double plot_x_spread = 0.06, plot_linear_width = 0.1;
inline constexpr double plot_ratio_min = 0, plot_ratio_max = 0; // Both zero: automatic.

// External prefix positions are one-based. Defaults follow the DES catalog.
inline constexpr auto external_columns = [] {
    bcg::ExternalColumns cols{};
    cols.ra = 5;
    cols.dec = 6;
    cols.zp = 17;
    cols.mag = 11;
    cols.flags_footprint = 1;
    cols.flags_foreground = 2;
    cols.flags_gold = 3;
    cols.ext_mash = 4;
    return cols;
}();

// ==========================================
// Function: Default source cuts initializer
// Method: Preserve the DES thresholds while exposing the full set in one file.
// ==========================================
inline constexpr auto cuts = [] {
    bcg::Cuts c{};
    c.snrfcut = 0;
    c.snrlow = 0;
    c.snrhigh = 0;
    c.starcut = 16;
    c.chi2_thresh = 0.01;
    c.flagcut = 0;
    c.zplow = 0;
    c.zphigh = 5;
    c.gf_lim = 0.0015;
    c.imaxcut = 64;
    c.jmaxcut = 64;
    c.mag_low = 0;
    c.mag_high = 25;
    c.psf_chi2_mltp = 3;
    c.r_half_thresh = 0;
    c.ft_cut = -1;
    c.fg_cut = -10;
    c.gold_cut = -10;
    c.ext_cut = -4;
    c.star_bar_mltp = -3;
    c.expo_bar_mltp = -3;
    c.max_dup = 5;
    return c;
}();
} // namespace bcg_esd::defaults
