#pragma once
#include "bcg_calibration.hpp"
#include "bcg_shear_cuts.hpp"

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
inline constexpr double gf_lim = -9999.0;
inline constexpr double selfcal_fd_limit = 0.0015; // When the source FD cut is disabled.
inline constexpr bcg::ShearCalibration calibration{0.0, 0.0, 0.0, 0.0}; // m1, m2, c1, c2

// Plot paths default to the result stem plus .png and _FD.png.
inline constexpr char esd_plot[] = "", fd_plot[] = "";
inline constexpr char plot_title[] = "BCG-centered group lensing";
inline constexpr double plot_x_spread = 0.06, plot_linear_width = 0.1;
inline constexpr double plot_ratio_min = 0, plot_ratio_max = 0; // Both zero: automatic.

// Positions are one-based. CatalogCombiner appends the exposure chi2.
inline constexpr auto external_columns = [] {
    bcg::ExternalColumns cols{};
    cols.ra = 5;
    cols.dec = 6;
    cols.zp = 17;
    cols.mag = 11;
    cols.ext_mash = 4;
    cols.expo_num = 0; // Set explicitly for the per-exposure bar.
    cols.expo_chi2 = -1; // Last column.
    return cols;
}();

// ==========================================
// Function: Default source cuts initializer
// Method: Start with every catalog cut disabled; the INI can supply measurement thresholds.
// ==========================================
inline constexpr auto cuts = [] {
    bcg::Cuts c{};
    c.snrfcut = -9999.0;
    c.snrlow = -9999.0;
    c.snrhigh = -9999.0;
    c.starcut = -9999.0;
    c.chi2_thresh = -9999.0;
    c.flagcut = -9999.0;
    c.zplow = -9999.0;
    c.zphigh = -9999.0;
    c.gf_lim = gf_lim;
    c.imaxcut = -9999.0;
    c.jmaxcut = -9999.0;
    c.mag_low = -9999.0;
    c.mag_high = -9999.0;
    c.psf_chi2_mltp = -9999.0;
    c.r_half_thresh = -9999.0;
    c.ext_cut = -9999.0;
    c.star_bar_mltp = -9999.0;
    c.expo_bar_mltp = -9999.0;
    c.max_dup = -9999.0;
    return c;
}();
} // namespace bcg_esd::defaults
