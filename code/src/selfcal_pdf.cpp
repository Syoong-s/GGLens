#include "selfcalibration.hpp"
#include "bcg_pdf_sym.hpp"
#include "bcg_quantiles.hpp"

namespace bcg_esd {
namespace {
// ==========================================
// Function: fd_bin
// Method: Match cochko's last-edge overflow rule before scanning refined edges.
// ==========================================
int fd_bin(double value, const std::vector<double> &edges) {
    if (value > edges.back()) return static_cast<int>(edges.size());
    for (int k = 0; k < static_cast<int>(edges.size()); ++k)
        if (value <= edges[k]) return k;
    return static_cast<int>(edges.size());
}

// ==========================================
// Function: occupancy
// Method: Stream the same global occupancy residual used by the original FD program.
// ==========================================
double occupancy(const bcg::PairStore &store, const std::vector<double> &edges,
                 int bin, std::int64_t total, MPI_Comm comm) {
    std::int64_t local = 0, global = 0;
    store.scan([&](const bcg::PairRecord &p) { if (fd_bin(std::abs(p.gt), edges) == bin) ++local; });
    MPI_Allreduce(&local, &global, 1, MPI_INT64_T, MPI_SUM, comm);
    return static_cast<double>(global) / total - 1.0 / (edges.size() + 1.0);
}

// ==========================================
// Function: equal_edges
// Method: Use bounded exact local sorting, original rank averaging and sparse global fallback.
// ==========================================
std::vector<double> equal_edges(const bcg::PairStore &store, int nbin,
                               std::int64_t total, const bcg::ScratchSpace &scratch,
                               std::size_t sort_bytes, MPI_Comm comm) {
    int rank = 0, ranks = 1;
    MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
    std::vector<double> local(nbin, 0), sums(nbin, 0), edges(nbin, 0);
    int efficient = (rank != 0 || ranks == 1) && store.count >= nbin + 1 ? 1 : 0;
    if (efficient) {
        const auto exact = bcg::prepare_quantiles(store, {0}, 0, nbin, 1, scratch, sort_bytes);
        local = exact.values[0];
    }
    MPI_Reduce(local.data(), sums.data(), nbin, MPI_DOUBLE, MPI_SUM, 0, comm);
    int n_eff = 0;
    MPI_Allreduce(&efficient, &n_eff, 1, MPI_INT, MPI_SUM, comm);
    if (n_eff == 0) {
        // Each normal worker has fewer than nbin+1 rows in this fallback.
        // An unusual nonempty manager is guarded by the same explicit memory quota.
        if (store.count > std::numeric_limits<int>::max() ||
            static_cast<std::uint64_t>(store.count) > sort_bytes / sizeof(double))
            throw std::runtime_error("Sparse self-calibration quantiles exceed rank memory budget");
        std::vector<double> values;
        values.reserve(static_cast<std::size_t>(store.count));
        store.scan([&](const bcg::PairRecord &p) { values.push_back(std::abs(p.gt)); });
        const int n = static_cast<int>(values.size());
        std::vector<int> counts(ranks), offsets(ranks);
        MPI_Allgather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
        std::int64_t gathered_n = 0;
        for (int i = 0; i < ranks; ++i) {
            if (gathered_n > std::numeric_limits<int>::max() - counts[i])
                throw std::runtime_error("Sparse FD gather exceeds MPI count limit");
            offsets[i] = static_cast<int>(gathered_n); gathered_n += counts[i];
        }
        if (static_cast<std::uint64_t>(gathered_n) > sort_bytes / sizeof(double))
            throw std::runtime_error("Sparse self-calibration gather exceeds memory budget");
        std::vector<double> all(rank == 0 ? static_cast<std::size_t>(gathered_n) : 0);
        MPI_Gatherv(values.data(), n, MPI_DOUBLE, all.data(), counts.data(), offsets.data(), MPI_DOUBLE, 0, comm);
        if (!rank) {
            std::sort(all.begin(), all.end());
            const auto dn = total / (nbin + 1);
            for (int i = 0; i < nbin; ++i) edges[i] = all.at(static_cast<std::size_t>(dn * (i + 1) - 1));
        }
    } else if (!rank) {
        for (int i = 0; i < nbin; ++i) edges[i] = sums[i] / n_eff;
    }
    MPI_Bcast(edges.data(), nbin, MPI_DOUBLE, 0, comm);
    for (int i = 0; i < nbin; ++i) {
        const double vc = occupancy(store, edges, i, total, comm);
        if (std::abs(vc) < 0.03) continue;
        double lo, hi;
        if (vc > 0) {
            lo = i == 0 ? (i + 1 < nbin ? 1.5*edges[i] - 0.5*edges[i+1] : 0.5*edges[i])
                        : 0.5 * (edges[i] + edges[i-1]);
            hi = edges[i];
        } else {
            lo = edges[i];
            hi = i + 1 >= nbin ? (i > 0 ? 1.5*edges[i] - 0.5*edges[i-1] : 1.5*edges[i])
                               : 0.5 * (edges[i] + edges[i+1]);
        }
        for (int iter = 0; iter < 20; ++iter) {
            edges[i] = lo;
            const double v1 = occupancy(store, edges, i, total, comm);
            if (std::abs(v1) < 0.03) break;
            edges[i] = hi;
            const double v2 = occupancy(store, edges, i, total, comm);
            if (std::abs(v2) < 0.03) break;
            if (v1*v2 < 0) {
                if (std::abs(v1) > std::abs(v2)) lo = 0.75*lo + 0.25*hi;
                else hi = 0.75*hi + 0.25*lo;
            } else if (v1 > 0) {
                const double old = lo; lo = 2.0*lo - hi; hi = 1.3*old - 0.3*hi;
            } else {
                const double old = hi; hi = 2.0*hi - lo; lo = 1.3*old - 0.3*lo;
            }
        }
    }
    return edges;
}

// ==========================================
// Function: chi2_grid
// Method: Scan bounded records for all trials, counting zero residuals in FD occupancy.
// ==========================================
std::vector<double> chi2_grid(const bcg::PairStore &store, const std::vector<double> &centers,
                              const std::vector<double> &edges, MPI_Comm comm) {
    const std::size_t bins = edges.size() + 1, span = 2*bins;
    std::vector<std::int64_t> local(centers.size()*span, 0), global(local.size());
    store.scan([&](const bcg::PairRecord &p) {
        for (std::size_t t = 0; t < centers.size(); ++t) {
            const double residual = p.gt - centers[t]*p.det;
            const auto cell = t*span + 2*static_cast<std::size_t>(fd_bin(std::abs(residual), edges));
            ++local[cell];
            if (residual > 0) ++local[cell+1];
            else if (residual < 0) --local[cell+1];
        }
    });
    MPI_Allreduce(local.data(), global.data(), static_cast<int>(local.size()), MPI_INT64_T, MPI_SUM, comm);
    std::vector<double> result(centers.size(), 0);
    for (std::size_t t = 0; t < centers.size(); ++t) for (std::size_t b = 0; b < bins; ++b) {
        const auto cell = t*span + 2*b;
        if (global[cell] > 0) {
            const double diff = static_cast<double>(global[cell+1]);
            result[t] += diff*diff / (2.0*static_cast<double>(global[cell]));
        }
    }
    return result;
}
} // namespace

// ==========================================
// Function: selfcal_fd_pdf
// Method: Match cochko's zero-centered search, threshold-20 shrinking and 200-point fit.
// ==========================================
FdPdfFit selfcal_fd_pdf(const bcg::PairStore &store, int magnitude_bins, int min_count,
                       double search_width, const bcg::ScratchSpace &scratch,
                       std::size_t sort_bytes, MPI_Comm comm) {
    const std::int64_t local_n = store.count;
    std::int64_t total = 0;
    MPI_Allreduce(&local_n, &total, 1, MPI_INT64_T, MPI_SUM, comm);
    if (total < std::max(min_count, magnitude_bins+1)) return {};
    const auto edges = equal_edges(store, magnitude_bins-1, total, scratch, sort_bytes, comm);
    double left = -search_width, right = search_width;
    bool converged = false;
    for (int iter = 0; iter < 200; ++iter) {
        const double mid = 0.5*(left+right);
        const auto chi2 = chi2_grid(store, {0.5*(left+mid), 0.5*(right+mid), mid}, edges, comm);
        bool changed = false;
        if (chi2[0] > chi2[2]+20.0) { left = 0.5*(left+mid); changed = true; }
        if (chi2[1] > chi2[2]+20.0) { right = 0.5*(right+mid); changed = true; }
        if (!changed) { converged = true; break; }
        if (!(right > left)) return {};
    }
    if (!converged) return {};
    std::vector<double> grid(200);
    for (int i = 0; i < 200; ++i) grid[i] = left+(right-left)*i/199;
    const auto values = chi2_grid(store, grid, edges, comm);
    const auto coeff = bcg::solve_quadratic(grid, values);
    if (!(coeff[0] > 0) || !std::isfinite(coeff[0])) return {};
    const double value = -coeff[1]/(2.0*coeff[0]);
    const double sigma = 1.0/std::sqrt(2.0*coeff[0]);
    if (!std::isfinite(value) || !std::isfinite(sigma) || sigma <= 0 || value <= left || value >= right) return {};
    return {value, sigma, true};
}
} // namespace bcg_esd
