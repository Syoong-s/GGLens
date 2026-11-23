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
    out.low = low;
    out.high = high;
    out.tangential_replicas.assign(stack.regions, out.measured);
    out.cross_replicas.assign(stack.regions, out.measured);
    MPI_Allreduce(&pairs.count, &out.count, 1, MPI_INT64_T, MPI_SUM, comm);
    if (out.count == 0) return out;

    double global_logr = 0;
    MPI_Allreduce(&pairs.log_radius_sum, &global_logr, 1, MPI_DOUBLE, MPI_SUM, comm);
    out.effective = std::exp(global_logr / out.count);
    if (out.count < cfg.pdf_bins + 2LL) return out;

    pairs.load_if_small(cfg.memory_bytes() / 2);
    const auto measure = [&](int component, double &full, std::vector<double> &replicas) {
        const auto pdf = bcg::prepare_replica_pdf(pairs, stack.labels, stack.regions,
                                                  component, cfg.pdf_bins, comm);
        const auto estimates = bcg::fit_replica_pdf(pairs, stack.labels, component, pdf, comm);
        full = estimates[0];
        replicas.assign(estimates.begin() + 1, estimates.end());
        if (component == 0) {
            std::vector<std::int64_t> omitted;
            for (int r = 0; r < stack.regions; ++r) omitted.push_back(out.count - pdf.kept[r + 1]);
            int rank = 0;
            MPI_Comm_rank(comm, &rank);
            if (!rank) {
                const auto limits = std::minmax_element(omitted.begin(), omitted.end());
                std::cout << "  Pairs removed per region: " << *limits.first << " to "
                          << *limits.second << std::endl;
            }
        }
    };
    measure(0, out.measured, out.tangential_replicas);
    measure(1, out.cross, out.cross_replicas);
    return out;
}
} // namespace bcg_esd
