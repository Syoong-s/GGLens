#include "bcg_mpi_rows.hpp"
#include "esd_fit_backend.hpp"
#include <iostream>

namespace bcg_esd {
namespace {
// ==========================================
// Function: fit_input
// Method: Extract the supported radial submatrix in the same order as the observation vector.
// ==========================================
FitInput fit_input(const std::vector<Result> &results, const Covariance &covariance) {
    FitInput input;
    for (const int i : covariance.fit_bins) {
        input.radius.push_back(results[i].effective);
        input.measured.push_back(results[i].measured);
        for (const int j : covariance.fit_bins)
            input.covariance.push_back(covariance.matrix[i*results.size() + j]);
    }
    return input;
}

// ==========================================
// Function: gather_chain
// Method: Gather ordered ensemble blocks only to the rank writing final products.
// ==========================================
std::vector<double> gather_chain(const std::vector<double> &local, int ensembles,
                                 int stride, MPI_Comm comm) {
    int rank = 0, ranks = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    std::vector<int> counts(ranks), offsets(ranks);
    for (int r = 0; r < ranks; ++r) {
        const auto [begin, end] = bcg::partition(ensembles, r, ranks);
        offsets[r] = static_cast<int>(begin * stride);
        counts[r] = static_cast<int>((end - begin) * stride);
    }
    if (local.size() != static_cast<std::size_t>(counts[rank]))
        throw std::runtime_error("Invalid local emcee chain shape");
    std::vector<double> combined(rank == 0 ? static_cast<std::size_t>(ensembles) * stride : 0);
    MPI_Gatherv(local.data(), counts[rank], MPI_DOUBLE, combined.data(),
                counts.data(), offsets.data(), MPI_DOUBLE, 0, comm);
    return combined;
}
} // namespace

// ==========================================
// Function: fit_esd
// Method: Share one CAMB run, sample independent ensembles on MPI ranks, and combine retained draws.
// ==========================================
void fit_esd(const std::vector<Result> &results, const Covariance &covariance,
             const LensStack &stack, const Config &cfg, MPI_Comm comm) {
    const auto fit_bins = covariance.fit_bins.size();
    if (fit_bins <= 5)
        throw std::runtime_error("Five-parameter ESD fit needs at least six complete jackknife bins");
    if (static_cast<std::size_t>(stack.regions) <= fit_bins)
        throw std::runtime_error("A full-rank covariance needs more jackknife regions than fitted radial bins");

    int rank = 0, ranks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    const auto [begin, end] = bcg::partition(cfg.fit_ensembles, rank, ranks);
    std::unique_ptr<FitBackend> backend;
    if (begin < end) backend = std::make_unique<FitBackend>(cfg);
    std::vector<double> corr;
    if (!rank) {
        std::cout << "Computing nonlinear CAMB/Halofit correlation at z=" << stack.redshift << std::endl;
        corr = backend->correlation(cfg, stack);
    }
    int size = static_cast<int>(corr.size());
    MPI_Bcast(&size, 1, MPI_INT, 0, comm);
    corr.resize(size);
    MPI_Bcast(corr.data(), size, MPI_DOUBLE, 0, comm);

    FitChain local;
    const int stride = 6 * cfg.fit_walkers * (cfg.fit_steps - cfg.fit_burn);
    if (backend) {
        backend->setup(fit_input(results, covariance), corr, cfg, stack);
        local.samples.reserve((end - begin) * stride);
        local.diagnostics.reserve((end - begin) * 6);
        for (auto ensemble = begin; ensemble < end; ++ensemble) {
            auto chain = backend->sample(cfg, static_cast<int>(ensemble));
            local.samples.insert(local.samples.end(), chain.samples.begin(), chain.samples.end());
            local.diagnostics.insert(local.diagnostics.end(), chain.diagnostics.begin(), chain.diagnostics.end());
            std::cout << "Rank " << rank << " sampled ensemble " << ensemble + 1
                      << '/' << cfg.fit_ensembles << std::endl;
        }
    }
    FitChain combined{gather_chain(local.samples, cfg.fit_ensembles, stride, comm),
                      gather_chain(local.diagnostics, cfg.fit_ensembles, 6, comm)};
    if (!rank) backend->finish(combined, cfg, stack, covariance.fit_bins);
}
} // namespace bcg_esd
