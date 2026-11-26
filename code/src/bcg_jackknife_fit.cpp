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

// Each component retains its own search range and minima across shared scans.
struct ReplicaSearch {
    const ReplicaPdf &pdf;
    std::vector<double> value, center, grid;
    std::vector<int> minima;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    bool finished;
    static constexpr int coarse_size = 41;

    explicit ReplicaSearch(const ReplicaPdf &pdf) : pdf(pdf),
        value(pdf.edges.size(), std::numeric_limits<double>::quiet_NaN()),
        center(pdf.edges.size(), 0), minima(pdf.edges.size(), -1),
        finished(pdf.edges[0].empty()) {
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
    }

    void coarse(const ReplicaCosts &chi2, int expansion) {
        bool low = false, high = false;
        for (std::size_t s = 0; s < pdf.edges.size(); ++s) {
            if (pdf.edges[s].empty()) continue;
            minima[s] = minimum(chi2[s], grid, center[s]);
            low |= minima[s] <= 1;
            high |= minima[s] >= coarse_size - 2;
        }
        if ((!low && !high) || expansion == 5) {
            finished = true;
            return;
        }
        const double width = hi - lo;
        if (low) lo -= width;
        if (high) hi += width;
    }

    void refine_grid() {
        const auto limits = std::minmax_element(minima.begin(), minima.end());
        int first = *limits.second, last = *limits.second;
        for (const auto m : minima) if (m >= 0) first = std::min(first, m);
        grid = trial_grid(grid[std::max(0, first - 3)], grid[std::min(coarse_size - 1, last + 3)], fine_trial_points);
    }

    void finish(const ReplicaCosts &chi2) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
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
    }
};
} // namespace

std::array<std::vector<double>, 2> fit_replica_pdf(const PairStore &pairs,
                                                const std::vector<int> &labels,
                                                const ComponentPdfs &pdf, MPI_Comm comm) {
    std::array<ReplicaSearch, 2> search{ReplicaSearch(pdf[0]), ReplicaSearch(pdf[1])};
    for (int expansion = 0; expansion <= 5; ++expansion) {
        TrialGrids grids{};
        for (int c = 0; c < 2; ++c) {
            auto &state = search[c];
            if (state.finished) continue;
            state.grid = trial_grid(state.lo, state.hi, ReplicaSearch::coarse_size);
            grids[c] = &state.grid;
        }
        if (!grids[0] && !grids[1]) break;
        const auto costs = replica_chi2(pairs, labels, pdf, grids, comm);
        for (int c = 0; c < 2; ++c)
            if (grids[c]) search[c].coarse(costs[c], expansion);
    }
    TrialGrids grids{};
    for (int c = 0; c < 2; ++c) {
        if (pdf[c].edges[0].empty()) continue;
        search[c].refine_grid();
        grids[c] = &search[c].grid;
    }
    const auto costs = replica_chi2(pairs, labels, pdf, grids, comm);
    for (int c = 0; c < 2; ++c)
        if (grids[c]) search[c].finish(costs[c]);
    return {std::move(search[0].value), std::move(search[1].value)};
}
} // namespace bcg
