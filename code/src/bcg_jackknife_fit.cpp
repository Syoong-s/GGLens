#include "bcg_jackknife.hpp"
#include "bcg_pdf_sym.hpp"

namespace bcg {
namespace {
// ==========================================
// Function: trial_grid
// Method: Evaluate a common ordered grid for the full sample and every deletion.
// ==========================================
std::vector<double> trial_grid(double lo, double hi, int size) {
    std::vector<double> grid(size);
    for (int i = 0; i < size; ++i) grid[i] = lo + (hi - lo) * i / (size - 1);
    return grid;
}

// ==========================================
// Function: minimum
// Method: Resolve tied minima towards the data-derived linear center.
// ==========================================
int minimum(const std::vector<double> &chi2, const std::vector<double> &grid, double center) {
    int best = 0;
    for (int i = 1; i < static_cast<int>(chi2.size()); ++i)
        if (chi2[i] < chi2[best] || (chi2[i] == chi2[best] &&
            std::abs(grid[i] - center) < std::abs(grid[best] - center))) best = i;
    return best;
}
} // namespace

// ==========================================
// Function: fit_replica_pdf
// Method: Bracket the sign-test minima and refine them locally without a halo model prior.
// ==========================================
std::vector<double> fit_replica_pdf(const PairStore &pairs, const std::vector<int> &labels,
                                    int component, const ReplicaPdf &pdf, MPI_Comm comm) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> value(pdf.edges.size(), nan), center(pdf.edges.size(), 0);
    if (pdf.edges[0].empty()) return value;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (std::size_t s = 0; s < pdf.edges.size(); ++s) {
        if (pdf.edges[s].empty()) continue;
        const auto [dd, dx, xx] = pdf.moments[s];
        center[s] = dx / dd;
        const double scatter = std::sqrt(std::max(0.0, xx - dx*center[s]) /
                                         (dd * (pdf.kept[s] - 1)));
        const double width = std::max({7 * scatter, 0.02 * std::abs(center[s]), 1e-6});
        lo = std::min(lo, center[s] - width);
        hi = std::max(hi, center[s] + width);
    }
    std::vector<double> grid;
    std::vector<int> minima(pdf.edges.size(), -1);
    constexpr int coarse_size = 41;
    for (int expansion = 0; expansion <= 5; ++expansion) {
        grid = trial_grid(lo, hi, coarse_size);
        const auto chi2 = replica_chi2(pairs, labels, component, pdf.edges, grid, comm);
        bool low = false, high = false;
        for (std::size_t s = 0; s < pdf.edges.size(); ++s) {
            if (pdf.edges[s].empty()) continue;
            minima[s] = minimum(chi2[s], grid, center[s]);
            low |= minima[s] <= 1;
            high |= minima[s] >= coarse_size - 2;
        }
        if ((!low && !high) || expansion == 5) break;
        const double width = hi - lo;
        if (low) lo -= width;
        if (high) hi += width;
    }
    const auto limits = std::minmax_element(minima.begin(), minima.end());
    int first = *limits.second, last = *limits.second;
    for (const auto m : minima) if (m >= 0) first = std::min(first, m);
    grid = trial_grid(grid[std::max(0, first - 3)], grid[std::min(coarse_size - 1, last + 3)], 61);
    const auto chi2 = replica_chi2(pairs, labels, component, pdf.edges, grid, comm);
    for (std::size_t s = 0; s < pdf.edges.size(); ++s) {
        if (minima[s] < 0) continue;
        const int best = minimum(chi2[s], grid, center[s]);
        if (best == 0 || best + 1 == static_cast<int>(grid.size())) continue;
        const int left = std::max(0, best - 4);
        const int right = std::min(static_cast<int>(grid.size()), best + 5);
        const auto coeff = solve_quadratic({grid.begin() + left, grid.begin() + right},
                                           {chi2[s].begin() + left, chi2[s].begin() + right});
        const double location = coeff[0] > 0 ? -coeff[1] / (2 * coeff[0]) : nan;
        value[s] = location >= grid[left] && location <= grid[right - 1] ? location : grid[best];
    }
    return value;
}
} // namespace bcg
