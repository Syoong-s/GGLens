#include "bcg_esd.hpp"
#include "bcg_jackknife.hpp"
#include <iostream>

namespace bcg_esd {
// ==========================================
// Function: measure_bin
// Method: Measure full and delete-one tangential and cross PDF_SYM estimates in one radial bin.
// ==========================================
Result measure_bin(bcg::PairStore &pairs, double low, double high, const LensStack &stack,
                   const Config &cfg, MPI_Comm comm) {
    Result out;
    const int regions = stack.sample.regions;
    out.low = low;
    out.high = high;
    out.tangential_replicas.assign(regions, out.measured);
    out.cross_replicas.assign(regions, out.measured);
    MPI_Allreduce(&pairs.count, &out.count, 1, MPI_INT64_T, MPI_SUM, comm);
    if (out.count == 0) return out;

    double global_logr = 0;
    MPI_Allreduce(&pairs.log_radius_sum, &global_logr, 1, MPI_DOUBLE, MPI_SUM, comm);
    out.effective = std::exp(global_logr / out.count);
    if (out.count < cfg.pdf_bins + 2LL) return out;

    // Both components share the histogram and cache allowance. The decoder
    // separately keeps only one source chunk and one read buffer.
    const auto budget = cfg.memory_bytes();
    const auto histogram_bytes = 2 * (static_cast<std::uint64_t>(regions) + 1) *
        bcg::pdf_histogram_bytes_per_slot(cfg.pdf_bins);
    if (histogram_bytes > budget / 2)
        throw std::runtime_error("Jackknife histograms exceed run.memory-mb budget");
    const auto cache_bytes = budget / 2 - static_cast<std::size_t>(histogram_bytes);
    pairs.prepare_scan(pairs.indexed() ? std::min(cache_bytes, budget / 8) : cache_bytes);
    const auto pdf = bcg::prepare_replica_pdf(pairs, stack.labels, regions, cfg.pdf_bins, comm);
    const auto estimates = bcg::fit_replica_pdf(pairs, stack.labels, pdf, comm);
    out.measured = estimates[0][0];
    out.cross = estimates[1][0];
    out.tangential_replicas.assign(estimates[0].begin() + 1, estimates[0].end());
    out.cross_replicas.assign(estimates[1].begin() + 1, estimates[1].end());
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (!rank) {
        const auto limits = std::minmax_element(pdf[0].kept.begin() + 1, pdf[0].kept.end());
        std::cout << "  Pairs removed per region: " << out.count - *limits.second
                  << " to " << out.count - *limits.first << std::endl;
    }
    return out;
}
} // namespace bcg_esd
