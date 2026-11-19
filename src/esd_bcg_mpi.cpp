#include "bcg_esd.hpp"
#include "selfcalibration.hpp"
#include <iostream>

namespace {
// ==========================================
// Function: run
// Method: Run DR9 matching, source pairing, jackknife measurement, MCMC, and the final figure.
// ==========================================
void run(bcg_esd::Config &cfg, MPI_Comm comm) {
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
    if (!rank) {
        std::cout << "Lens center: " << bcg::center_name(cfg.center) << std::endl;
        if (cfg.self_calibrate) std::cout << "Self-calibration enabled: Cartesian all-R union; manual m/c will be overridden" << std::endl;
        else std::cout << "Shear calibration: " << cfg.calibration.description() << std::endl;
    }
    std::vector<std::string> files;
    if (!rank) files = bcg::read_shear_list(cfg.shear_list);
    bcg::broadcast_paths(files, comm);
    const auto layout = bcg::read_shear_layout(files.front(), cfg.external_columns, cfg.cuts, comm);
    if (!rank) std::cout << "First shear catalog columns: " << layout.destinations.size() << '\n'
                         << "Restore catalog FD: " << cfg.fd_corrected << '\n'
                         << "Shear files: " << files.size() << ", ranks: " << size
                         << ", buffer budget per rank: " << cfg.memory_mib << " MiB\n"
                         << "Scratch root: " << cfg.scratch_dir << std::endl;
    const auto lenses = bcg_esd::load_lenses(cfg, comm);
    bcg::ScratchSpace scratch(cfg.scratch_dir, comm);
    bcg_esd::SelfCalibrationDiagnostics diagnostics;
    auto bins = bcg_esd::build_pair_spools(lenses, files, layout, cfg, diagnostics, scratch, comm);
    const auto stack = bcg_esd::make_lens_stack(lenses, bins, cfg, comm);
    const auto edges = bcg_esd::radial_edges(cfg);
    std::vector<bcg_esd::Result> results; results.reserve(cfg.n_r_bins);
    for (int i = 0; i < cfg.n_r_bins; ++i) {
        if (!rank) std::cout << "Measuring bin " << i + 1 << '/' << cfg.n_r_bins
                             << " R=[" << edges[i] << ',' << edges[i + 1] << ")" << std::endl;
        auto result = bcg_esd::measure_bin(*bins[i], edges[i], edges[i + 1], stack, cfg, comm);
        bins[i]->release(); bins[i].reset();
        if (!rank) {
            std::cout << "bin " << i + 1 << '/' << cfg.n_r_bins << " pairs=" << result.count
                      << " R_eff=" << result.effective << " ESD=" << result.measured << std::endl;
        }
        results.push_back(std::move(result));
    }
    scratch.cleanup(comm);
    const auto covariance = bcg_esd::jackknife_covariance(results, stack.regions);
    if (!rank) {
        bcg_esd::write_results(results, covariance, stack, cfg);
        std::cout << "Saved ESD and jackknife covariance: " << cfg.output << std::endl;
    }
    bcg_esd::fit_esd(results, covariance, stack, cfg, comm);
    if (!rank) std::cout << "Fit figure: " << cfg.fit_plot << std::endl;
}
} // namespace

// ==========================================
// Function: main
// Method: Accept an optional INI path, run the measurement with MPI, and abort on error.
// ==========================================
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try {
        if (argc > 2)
            throw std::runtime_error(std::string("Usage: ") + argv[0] + " [ini-path]");
        const char *ini_path = argc == 2 ? argv[1] : nullptr;
        bcg_esd::Config cfg = bcg_esd::load_config(ini_path, MPI_COMM_WORLD);
        run(cfg, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Rank " << rank << " error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
