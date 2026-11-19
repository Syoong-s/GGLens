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
    bcg::Center center = defaults::center;
    double zmin = defaults::zmin, zmax = defaults::zmax;
    double mmin = defaults::mmin, mmax = defaults::mmax;
    double rmin = defaults::rmin, rmax = defaults::rmax, dzmin = defaults::dzmin;
    double omega_m = defaults::omega_m, h = defaults::h, sigma8 = defaults::sigma8;
    double omega_b_h2 = defaults::omega_b_h2, spectral_index = defaults::spectral_index;
    double neutrino_mass = defaults::neutrino_mass;
    int n_r_bins = defaults::n_r_bins, n_jack = defaults::n_jack, pdf_bins = defaults::pdf_bins;
    int memory_mib = defaults::memory_mib;
    bool fd_corrected = defaults::fd_corrected, self_calibrate = defaults::self_calibrate;
    int selfcal_fd_bins = defaults::selfcal_fd_bins, selfcal_pdf_bins = defaults::selfcal_pdf_bins;
    int selfcal_min_count = defaults::selfcal_min_count;
    double selfcal_search_width = defaults::selfcal_search_width;
    int fit_walkers = defaults::fit_walkers, fit_steps = defaults::fit_steps, fit_burn = defaults::fit_burn;
    int fit_ensembles = defaults::fit_ensembles, fit_seed = defaults::fit_seed;
    int fit_profile_draws = defaults::fit_profile_draws;
    std::string fit_python_path = defaults::fit_python_path;
    std::string fit_plot = defaults::fit_plot, plot_title = defaults::plot_title;
    std::string covariance_output, replica_output, summary_output, chain_output;
    std::string profile_output, profile_draws_output;
    bcg::ExternalColumns external_columns = defaults::external_columns;
    bcg::Cuts cuts{defaults::cuts, defaults::excluded_ccds};
    bcg::ShearCalibration calibration = defaults::calibration;
    std::int64_t self_calibration_sources = 0;
    // ==========================================
    // Function: fd_limit
    // Method: Keep the self-calibration grid defined when the catalog FD cut is disabled.
    // ==========================================
    double fd_limit() const {
        return bcg::cut_enabled(cuts.gf_lim) ? cuts.gf_lim : defaults::selfcal_fd_limit;
    }
    // ==========================================
    // Function: memory_bytes
    // Method: Convert the rank-local MiB budget to an unsigned byte count.
    // ==========================================
    std::size_t memory_bytes() const { return static_cast<std::size_t>(memory_mib) * 1024 * 1024; }
};
struct LensStack {
    std::vector<int> labels;
    std::size_t n_lenses = 0;
    int regions = 0;
    double stellar_mass = 0; // Mean in h^-1 M_sun.
    double redshift = 0;
};
struct Result {
    double low = 0, high = 0, effective = std::numeric_limits<double>::quiet_NaN();
    std::int64_t count = 0;
    double measured = std::numeric_limits<double>::quiet_NaN(), sigma = measured;
    double cross = measured, cross_sigma = measured;
    std::vector<double> tangential_replicas, cross_replicas;
};
struct Covariance {
    std::vector<double> matrix; // Row-major in radial-bin order.
    std::vector<int> fit_bins;
};
struct SelfCalibrationDiagnostics;
Config load_config(const char *ini_path, MPI_Comm comm);
std::vector<bcg::Lens> load_lenses(const Config &cfg, MPI_Comm comm);
LensStack make_lens_stack(const std::vector<bcg::Lens> &lenses,
                          const std::vector<std::unique_ptr<bcg::PairStore>> &bins,
                          const Config &cfg, MPI_Comm comm);
std::vector<double> radial_edges(const Config &cfg);
std::vector<std::unique_ptr<bcg::PairStore>> build_pair_spools(
    const std::vector<bcg::Lens> &lenses, const std::vector<std::string> &files,
    const bcg::ShearLayout &layout, Config &cfg, SelfCalibrationDiagnostics &diagnostics,
    const bcg::ScratchSpace &scratch, MPI_Comm comm);
Result measure_bin(bcg::PairStore &pairs, double low, double high, const LensStack &stack,
                   const Config &cfg, MPI_Comm comm);
Covariance jackknife_covariance(std::vector<Result> &results, int regions);
void write_results(const std::vector<Result> &results, const Covariance &covariance,
                   const LensStack &stack, const Config &cfg);
void fit_esd(const std::vector<Result> &results, const Covariance &covariance,
             const LensStack &stack, const Config &cfg, MPI_Comm comm);
} // namespace bcg_esd
