#include "bcg_jackknife.hpp"

namespace bcg {
namespace {
// ==========================================
// Function: add_trial_intervals
// Method: Add one residual to contiguous trial ranges with the same sign and magnitude bin.
// ==========================================
void add_trial_intervals(double value, double response, const std::vector<double> &trials,
                         const std::vector<double> &edges, std::int64_t *counts, std::size_t stride) {
    const auto residual = [&](std::size_t trial) { return value - trials[trial] * response; };
    const auto code = [&](std::size_t trial) {
        const double r = residual(trial);
        if (r == 0) return 0;
        const int bin = 1 + static_cast<int>(std::lower_bound(edges.begin(), edges.end(),
            std::abs(r)) - edges.begin());
        return r > 0 ? bin : -bin;
    };

    std::size_t begin = 0;
    int current = code(begin);
    const double last = residual(trials.size() - 1);
    while (begin < trials.size()) {
        const int bin = std::abs(current) - 1;
        const double lower = bin > 0 ? edges[bin - 1] : 0;
        const double upper = bin >= 0 && static_cast<std::size_t>(bin) < edges.size()
            ? edges[bin] : std::numeric_limits<double>::infinity();
        const auto in_bin = [&](double r) {
            if (current == 0) return r == 0;
            if (std::isnan(r)) return current == -1;
            const double signed_magnitude = current > 0 ? r : -r;
            return signed_magnitude > lower && signed_magnitude <= upper;
        };

        auto end = trials.size();
        if (!in_bin(last)) {
            // The signed bin is monotone along the grid. Testing its two
            // boundaries gives the same result as a fresh edge lookup.
            auto lo = begin + 1, hi = trials.size() - 1;
            while (lo < hi) {
                const auto mid = lo + (hi - lo) / 2;
                if (in_bin(residual(mid))) lo = mid + 1;
                else hi = mid;
            }
            end = lo;
        }
        if (current != 0) {
            const int signed_bin = 2 * bin + (current < 0);
            ++counts[begin*stride + signed_bin];
            if (end < trials.size()) --counts[end*stride + signed_bin];
        }
        begin = end;
        if (begin < trials.size()) current = code(begin);
    }
}

class TrialAccumulator {
    const PdfEdges &edges;
    const std::vector<double> &trials;
    std::size_t regions, bins, union_bins = 0, stride = 0;
    std::vector<double> union_edges;
    std::vector<std::vector<std::size_t>> cdf_index;
    std::vector<std::int64_t> counts;

public:
    TrialAccumulator(const PdfEdges &edges, const std::vector<double> &trials)
        : edges(edges), trials(trials), regions(edges.size() - 1), bins(edges[0].size() + 1) {
        // With near-identical delete-one edges, a pair still updates only two bins per trial.
        for (const auto &row : edges) union_edges.insert(union_edges.end(), row.begin(), row.end());
        std::sort(union_edges.begin(), union_edges.end());
        union_edges.erase(std::unique(union_edges.begin(), union_edges.end()), union_edges.end());
        union_bins = union_edges.size() + 1;
        cdf_index.resize(edges.size());
        for (std::size_t s = 0; s < edges.size(); ++s) {
            if (edges[s].empty()) continue;
            for (const double edge : edges[s])
                cdf_index[s].push_back(std::lower_bound(union_edges.begin(), union_edges.end(), edge)
                                       - union_edges.begin());
            cdf_index[s].push_back(union_bins - 1);
        }
        stride = 2 * (union_bins + regions * bins);
        if (stride > static_cast<std::size_t>(std::numeric_limits<int>::max()) / trials.size())
            throw std::runtime_error("Jackknife sign histogram exceeds MPI count limit");
        counts.resize(trials.size() * stride);
    }

    void append(std::size_t region, double value, double response) {
        const auto &own_edges = edges[region + 1];
        const auto own_offset = 2 * (union_bins + region * bins);
        add_trial_intervals(value, response, trials, union_edges, counts.data(), stride);
        if (!own_edges.empty())
            add_trial_intervals(value, response, trials, own_edges, counts.data() + own_offset, stride);
    }

    ReplicaCosts finish(MPI_Comm comm) {
        MPI_Allreduce(MPI_IN_PLACE, counts.data(), static_cast<int>(counts.size()), MPI_INT64_T, MPI_SUM, comm);
        for (std::size_t j = 1; j < trials.size(); ++j)
            for (std::size_t i = 0; i < stride; ++i) counts[j*stride + i] += counts[(j-1)*stride + i];
        std::vector<std::vector<double>> chi2(edges.size(), std::vector<double>(trials.size()));
        for (std::size_t j = 0; j < trials.size(); ++j) {
            auto *counts_at_trial = counts.data() + j * stride;
            for (std::size_t b = 1; b < union_bins; ++b) {
                counts_at_trial[2*b] += counts_at_trial[2*(b-1)];
                counts_at_trial[2*b+1] += counts_at_trial[2*(b-1)+1];
            }
            for (std::size_t s = 0; s < edges.size(); ++s) {
                if (edges[s].empty()) continue;
                const auto *omitted = s > 0 ? counts_at_trial + 2*(union_bins + (s-1)*bins) : nullptr;
                std::int64_t prev_pos = 0, prev_neg = 0;
                for (std::size_t b = 0; b < bins; ++b) {
                    const auto at = cdf_index[s][b];
                    const auto ending_pos = counts_at_trial[2*at], ending_neg = counts_at_trial[2*at+1];
                    const double pos = ending_pos - prev_pos - (omitted ? omitted[2*b] : 0);
                    const double neg = ending_neg - prev_neg - (omitted ? omitted[2*b+1] : 0);
                    if (pos + neg > 0) chi2[s][j] += (pos - neg) * (pos - neg) / (2 * (pos + neg));
                    prev_pos = ending_pos;
                    prev_neg = ending_neg;
                }
            }
        }
        return chi2;
    }
};
} // namespace

std::array<ReplicaCosts, 2> replica_chi2(const PairStore &pairs, const std::vector<int> &labels,
                                        const ComponentPdfs &pdf, const TrialGrids &grids, MPI_Comm comm) {
    std::array<std::unique_ptr<TrialAccumulator>, 2> accumulators;
    for (int c = 0; c < 2; ++c)
        if (grids[c]) accumulators[c] = std::make_unique<TrialAccumulator>(pdf[c].edges, *grids[c]);
    std::array<ReplicaCosts, 2> results;
    if (!grids[0] && !grids[1]) return results;
    pairs.scan([&](const PairRecord &p) {
        const auto region = static_cast<std::size_t>(labels[p.lens]);
        if (accumulators[0]) accumulators[0]->append(region, p.gt, p.det);
        if (accumulators[1]) accumulators[1]->append(region, p.gx, p.dex);
    });
    for (int c = 0; c < 2; ++c) {
        if (!accumulators[c]) continue;
        results[c] = accumulators[c]->finish(comm);
        accumulators[c].reset();
    }
    return results;
}
} // namespace bcg
