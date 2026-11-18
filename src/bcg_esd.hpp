#pragma once
#include "bcg_catalog.hpp"
#include "esd_defaults.hpp"
#include "bcg_pair_store.hpp"

namespace bcg_esd {
// ==========================================
// Structure: Config
// Method: Start from compiled defaults and retain the runtime overrides and effective calibration.
// ==========================================
struct Config {
    std::string group_list = defaults::group_list, shear_list = defaults::shear_list;
    std::string output = defaults::output, cache = defaults::cache;
    std::string scratch_dir = defaults::scratch_dir;
    double zmin = defaults::zmin, zmax = defaults::zmax;
    double mmin = defaults::mmin, mmax = defaults::mmax;
    double rmin = defaults::rmin, rmax = defaults::rmax, dzmin = defaults::dzmin;
    double omega_m = defaults::omega_m, h = defaults::h, nfw_c = defaults::nfw_c;
    double conc_min = defaults::conc_min, conc_max = defaults::conc_max;
    std::string conc_model = defaults::conc_model;
    int n_r_bins = defaults::n_r_bins, n_jack = defaults::n_jack, pdf_bins = defaults::pdf_bins;
    int memory_mib = defaults::memory_mib;
    bool fd_corrected = defaults::fd_corrected, self_calibrate = defaults::self_calibrate;
    int selfcal_fd_bins = defaults::selfcal_fd_bins, selfcal_pdf_bins = defaults::selfcal_pdf_bins;
    int selfcal_min_count = defaults::selfcal_min_count;
    double selfcal_search_width = defaults::selfcal_search_width;
    std::string esd_plot = defaults::esd_plot, fd_plot = defaults::fd_plot;
    std::string plot_title = defaults::plot_title;
    double plot_x_spread = defaults::plot_x_spread, plot_linear_width = defaults::plot_linear_width;
    double plot_ratio_min = defaults::plot_ratio_min, plot_ratio_max = defaults::plot_ratio_max;
    bcg::ExternalColumns external_columns = defaults::external_columns;
    bcg::Cuts cuts = defaults::cuts;
    bcg::ShearCalibration calibration = defaults::calibration;
    std::int64_t self_calibration_sources = 0;
    // ==========================================
    // Function: memory_bytes
    // Method: Convert the rank-local MiB budget to an unsigned byte count.
    // ==========================================
    std::size_t memory_bytes() const { return static_cast<std::size_t>(memory_mib) * 1024 * 1024; }
};
struct Result {
    double low = 0, high = 0, effective = 0;
    std::int64_t count = 0;
    double measured = 0, sigma = 0, cross = 0, cross_sigma = 0;
    double theory = 0, ratio = 0, residual = 0;
};
struct SelfCalibrationDiagnostics;
Config load_config(MPI_Comm comm);
std::vector<bcg::Lens> load_lenses(const Config &cfg, MPI_Comm comm);
std::vector<int> initial_regions(const std::vector<bcg::Lens> &lenses, int regions, MPI_Comm comm);
std::vector<double> radial_edges(const Config &cfg);
std::vector<std::unique_ptr<bcg::PairStore>> build_pair_spools(
    const std::vector<bcg::Lens> &lenses, const std::vector<int> &labels,
    const std::vector<std::string> &files, const bcg::ShearLayout &layout, Config &cfg,
    SelfCalibrationDiagnostics &diagnostics, const bcg::ScratchSpace &scratch, MPI_Comm comm);
Result measure_bin(bcg::PairStore &pairs, double low, double high,
                   const std::vector<bcg::Lens> &lenses, const std::vector<int> &labels,
                   const Config &cfg, const bcg::ScratchSpace &scratch, MPI_Comm comm);
void write_results(const std::vector<Result> &results, const Config &cfg);
} // namespace bcg_esd
