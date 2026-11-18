#include "bcg_pdf_stream.hpp"

namespace bcg {
namespace {
// ==========================================
// Function: source_accumulate
// Method: Count the original absolute-shear occupancy through the selected pair view.
// ==========================================
double source_accumulate(const PairView &view, const std::vector<double> &edges,
                         int bin, std::int64_t ntot_ref, MPI_Comm comm) {
    std::int64_t local = 0, global = 0;
    view.scan([&](double value, double) { if (f77_bin(std::abs(value), edges) == bin) ++local; });
    MPI_Allreduce(&local, &global, 1, MPI_INT64_T, MPI_SUM, comm);
    return ntot_ref > 0 ? static_cast<double>(global) / ntot_ref - 1.0 / (edges.size() + 1.0) : 0;
}

// ==========================================
// Function: equal_prob_edges
// Method: Preserve rank quantile averaging and original occupancy refinement using views.
// ==========================================
std::vector<double> equal_prob_edges(const PairView &view, const std::vector<double> &initial,
                                     int nbin, std::int64_t ntot_ref, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    std::vector<double> local(nbin, 0), sums(nbin, 0), edges(nbin, 0);
    int efficient = rank != 0 && view.size() >= static_cast<std::int64_t>(nbin) + 1 ? 1 : 0;
    if (efficient) local = initial;
    MPI_Reduce(local.data(), sums.data(), nbin, MPI_DOUBLE, MPI_SUM, 0, comm);
    int n_eff = 0; MPI_Reduce(&efficient, &n_eff, 1, MPI_INT, MPI_SUM, 0, comm);
    if (!rank) {
        if (n_eff > 0) for (int i = 0; i < nbin; ++i) edges[i] = sums[i] / n_eff;
        else if (view.size() >= static_cast<std::int64_t>(nbin) + 1) edges = initial;
        else for (int i = 0; i < nbin; ++i)
            edges[i] = 0.001 + (0.1 - 0.001) * i / std::max(nbin - 1, 1);
    }
    MPI_Bcast(edges.data(), nbin, MPI_DOUBLE, 0, comm);

    for (int i = 0; i < nbin; ++i) {
        double vc = source_accumulate(view, edges, i, ntot_ref, comm);
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
            const double v1 = source_accumulate(view, edges, i, ntot_ref, comm);
            if (std::abs(v1) < 0.03) break;
            edges[i] = hi;
            const double v2 = source_accumulate(view, edges, i, ntot_ref, comm);
            if (std::abs(v2) < 0.03) break;
            if (v1 * v2 < 0) {
                if (std::abs(v1) > std::abs(v2)) lo = 0.75*lo + 0.25*hi;
                else hi = 0.75*hi + 0.25*lo;
            } else if (v1 > 0) {
                const double old = lo;
                lo = 2.0*lo - hi;
                hi = 1.3*old - 0.3*hi;
            } else {
                const double old = hi;
                hi = 2.0*hi - lo;
                lo = 1.3*old - 0.3*lo;
            }
        }
    }
    return edges;
}

// ==========================================
// Function: chi2_grid
// Method: Scan a block once for all trial values and reduce only small sign histograms.
// ==========================================
std::vector<double> chi2_grid(const PairView &view, const std::vector<double> &centers,
                             const std::vector<double> &edges, MPI_Comm comm) {
    const int nb = static_cast<int>(edges.size()) + 1;
    const std::size_t span = 2 * static_cast<std::size_t>(nb);
    std::vector<std::int64_t> local(centers.size() * span, 0), global(local.size(), 0);
    view.scan([&](double x, double response) {
        for (std::size_t j = 0; j < centers.size(); ++j) {
            const double y = x - centers[j] * response;
            if (y == 0) continue;
            const int bin = f77_bin(std::abs(y), edges);
            ++local[j * span + (y > 0 ? bin : nb + bin)];
        }
    });
    MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()), MPI_INT64_T, MPI_SUM, comm);
    std::vector<double> result(centers.size(), 0);
    for (std::size_t j = 0; j < centers.size(); ++j) for (int b = 0; b < nb; ++b) {
        const double pos = static_cast<double>(global[j * span + b]);
        const double neg = static_cast<double>(global[j * span + nb + b]);
        if (pos + neg > 0) result[j] += (pos - neg) * (pos - neg) / (2 * (pos + neg));
    }
    return result;
}
} // namespace

// ==========================================
// Function: pdf_sym_stream
// Method: Equal-probability sign test, bracket search, and quadratic minimum.
// ==========================================
PdfResult pdf_sym_stream(const PairView &view, const std::vector<double> &initial, int nbin,
                                 MPI_Comm comm, double c0, double dc0,
                                 std::int64_t ntot_ref) {
    std::int64_t local_n = view.size(), total_n = 0;
    MPI_Allreduce(&local_n, &total_n, 1, MPI_INT64_T, MPI_SUM, comm);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (total_n < nbin + 2) return {nan, nan, nan};
    const auto edges = equal_prob_edges(view, initial, nbin, ntot_ref > 0 ? ntot_ref : total_n, comm);
    double local_sums[2] = {}, global_sums[2] = {};
    view.scan([&](double x, double response) {
        local_sums[0] += response * response;
        local_sums[1] += x * response;
    });
    MPI_Allreduce(local_sums, global_sums, 2, MPI_DOUBLE, MPI_SUM, comm);
    if (global_sums[0] <= 0) return {0, 0, chi2_grid(view, {0}, edges, comm)[0]};
    const double linear = global_sums[1] / global_sums[0];
    const double sigma_linear = 1.0 / std::sqrt(global_sums[0]);
    if (!std::isfinite(c0) || !std::isfinite(dc0)) {
        c0 = linear;
        dc0 = std::max(5*sigma_linear, std::abs(c0)*0.5 + 1e-10);
    }
    const double dc = dc0 > 0 ? dc0 : std::max(5*sigma_linear, std::abs(c0)*0.5 + 1e-10);
    double left = c0 - dc, right = c0 + dc;

    for (int expand = 0; expand <= 6; ++expand) {
        std::vector<double> grid(100);
        for (int i = 0; i < 100; ++i) grid[i] = left + (right - left)*i/99;
        const auto y = chi2_grid(view, grid, edges, comm);
        const int best = static_cast<int>(std::min_element(y.begin(), y.end()) - y.begin());
        if (best >= 3 && best <= 96) break;
        if (expand == 6) break;
        const double width = 2*(right - left);
        if (best < 3) left = right - width;
        else right = left + width;
    }

    for (int iter = 0; iter < 200; ++iter) {
        const double mid = 0.5*(left + right);
        const auto y = chi2_grid(view,
            {0.5*(left + mid), 0.5*(right + mid), mid}, edges, comm);
        bool changed = false;
        if (y[0] > y[2] + 20.0) { left = 0.5*(left + mid); changed = true; }
        if (y[1] > y[2] + 20.0) { right = 0.5*(right + mid); changed = true; }
        if (!changed) break;
    }

    std::vector<double> grid(200);
    for (int i = 0; i < 200; ++i) grid[i] = left + (right - left)*i/199;
    const auto y = chi2_grid(view, grid, edges, comm);
    const int best = static_cast<int>(std::min_element(y.begin(), y.end()) - y.begin());
    const auto coeff = solve_quadratic(grid, y);
    double value = grid[best], sigma = nan;
    if (coeff[0] > 0 && std::isfinite(coeff[0])) {
        const double candidate = -coeff[1] / (2*coeff[0]);
        if (std::isfinite(candidate) && candidate >= left && candidate <= right) {
            value = candidate;
            sigma = 1.0 / std::sqrt(2*coeff[0]);
        }
    }
    return {value, sigma, chi2_grid(view, {value}, edges, comm)[0]};
}

} // namespace bcg
