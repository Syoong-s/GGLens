#include "esd_plot.hpp"
#include <iostream>

namespace {
// ==========================================
// Function: run
// Method: Broadcast the input schema, measure ESD, and render the outputs from rank-zero results.
// ==========================================
void run(bcg_esd::Config &cfg, MPI_Comm comm) {
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
    if (!rank) {
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
    const auto labels = bcg_esd::initial_regions(lenses, cfg.n_jack, comm);
    bcg::ScratchSpace scratch(cfg.scratch_dir, comm);
    bcg_esd::SelfCalibrationDiagnostics diagnostics;
    auto bins = bcg_esd::build_pair_spools(lenses, labels, files, layout, cfg, diagnostics, scratch, comm);
    const auto edges = bcg_esd::radial_edges(cfg);
    std::vector<bcg_esd::Result> results; results.reserve(cfg.n_r_bins);
    for (int i = 0; i < cfg.n_r_bins; ++i) {
        if (!rank) std::cout << "Measuring bin " << i + 1 << '/' << cfg.n_r_bins
                             << " R=[" << edges[i] << ',' << edges[i + 1] << ")" << std::endl;
        const auto result = bcg_esd::measure_bin(*bins[i], edges[i], edges[i + 1],
                                                lenses, labels, cfg, scratch, comm);
        bins[i]->release(); bins[i].reset();
        if (!rank) {
            results.push_back(result);
            std::cout << "bin " << i + 1 << '/' << cfg.n_r_bins << " pairs=" << result.count
                      << " R_eff=" << result.effective << " ESD=" << result.measured << std::endl;
        }
    }
    scratch.cleanup(comm);
    if (!rank) {
        bcg_esd::write_results(results, cfg);
        std::cout << "Saved: " << cfg.output << std::endl;
        bcg_esd::render_plots(results, cfg, diagnostics);
        std::cout << "ESD plot: " << cfg.esd_plot << std::endl;
        if (cfg.self_calibrate) std::cout << "FD plot: " << cfg.fd_plot << std::endl;
    }
}
} // namespace

// ==========================================
// Function: main
// Method: Initialize MPI, run the complete BCG measurement, abort all ranks on error.
// ==========================================
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try {
        if (argc != 1)
            throw std::runtime_error("Configure in esd_defaults.hpp or esd.ini and launch without arguments");
        bcg_esd::Config cfg = bcg_esd::load_config(MPI_COMM_WORLD);
        run(cfg, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Rank " << rank << " error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
