#include "bcg_esd.hpp"
#include "esd_fit_backend.hpp"
#include <iostream>

namespace {
// ==========================================
// Function: measure
// Method: Measure the ESD and publish the data needed to fit it on a later run.
// ==========================================
void measure(bcg_esd::Config &cfg, MPI_Comm comm) {
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
    bcg_esd::Measurement measurement;
    auto bins = bcg_esd::build_pair_spools(lenses, files, layout, cfg, scratch,
                                         measurement.sources_after_cuts, comm);
    auto stack = bcg_esd::make_lens_stack(lenses, bins, cfg, comm);
    const auto edges = bcg_esd::radial_edges(cfg);
    if (!rank) measurement.bins.reserve(cfg.n_r_bins);
    for (int i = 0; i < cfg.n_r_bins; ++i) {
        if (!rank) std::cout << "Measuring bin " << i + 1 << '/' << cfg.n_r_bins
                             << " R=[" << edges[i] << ',' << edges[i + 1] << ")" << std::endl;
        auto result = bcg_esd::measure_bin(*bins[i], edges[i], edges[i + 1], stack, cfg, comm);
        bins[i]->release(); bins[i].reset();
        if (!rank) {
            std::cout << "bin " << i + 1 << '/' << cfg.n_r_bins << " pairs=" << result.count
                      << " R_eff=" << result.effective << " ESD=" << result.measured << std::endl;
            measurement.bins.push_back(std::move(result));
        }
    }
    scratch.cleanup(comm);
    if (rank) return;
    measurement.metadata = {cfg.center, cfg.cosmology, stack.sample, cfg.n_r_bins};
    measurement.groups_per_region = std::move(stack.groups_per_region);
    measurement.photometry_selected = stack.photometry_selected;
    measurement.requested_regions = cfg.n_jack;
    measurement.fd_restored = cfg.fd_corrected;
    measurement.calibration = cfg.calibration;
    if (cfg.self_calibrate) {
        measurement.calibration_sources = cfg.self_calibration_sources;
        measurement.fd_bins = cfg.selfcal_fd_bins;
        measurement.fd_limit = cfg.fd_limit();
    }
    bcg_esd::write_measurement(cfg.measurement_output, measurement);
    std::cout << "Saved ESD measurement: " << cfg.measurement_output << std::endl;
}

void run(bcg_esd::Config &cfg, MPI_Comm comm) {
    using bcg_esd::Stage;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (cfg.stage == Stage::pull) {
        if (!rank) bcg_esd::FitBackend(cfg).pull(cfg);
        return;
    }
    if (cfg.stage == Stage::all || cfg.stage == Stage::measure) measure(cfg, comm);
    if (cfg.stage == Stage::measure) return;
    bcg_esd::FitInput input;
    if (!rank) input = bcg_esd::load_fit_input(cfg);
    auto backend = bcg_esd::fit_esd(std::move(input), cfg, comm);
    if (!rank && cfg.stage == Stage::all) backend->pull(cfg);
}
} // namespace

// ==========================================
// Function: main
// Method: Accept an INI path and optional stage override and abort all ranks on error.
// ==========================================
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try {
        const std::string usage = std::string("Usage: ") + argv[0] + " [--stage all|measure|fit|pull] [ini-path]";
        const char *ini_path = nullptr;
        const char *stage = nullptr;
        for (int i = 1; i < argc; ++i) {
            const std::string option(argv[i]);
            if (option == "--help" || option == "-h") {
                if (!rank) std::cout << usage << std::endl;
                MPI_Finalize();
                return 0;
            }
            if (option == "--stage" && !stage && i + 1 < argc) {
                stage = argv[++i];
            } else if (!ini_path && !option.empty() && option.front() != '-') {
                ini_path = argv[i];
            } else {
                throw std::runtime_error(usage);
            }
        }
        bcg_esd::Config cfg = bcg_esd::load_config(ini_path, MPI_COMM_WORLD, stage);
        run(cfg, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Rank " << rank << " error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
