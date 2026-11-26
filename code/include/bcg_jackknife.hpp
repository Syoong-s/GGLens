#pragma once
#include "bcg_pair_store.hpp"

namespace bcg {
using PdfEdges = std::vector<std::vector<double>>;

inline constexpr int magnitude_histogram_cells = 2049;
inline constexpr int fine_trial_points = 61;

// A component needs one magnitude histogram or one trial histogram at a time.
inline std::uint64_t pdf_histogram_bytes_per_slot(int bins) {
    const auto trial_cells = (static_cast<std::uint64_t>(bins) + 1) * 4 * fine_trial_points;
    return std::max(trial_cells, static_cast<std::uint64_t>(magnitude_histogram_cells)) * sizeof(std::int64_t);
}

struct ReplicaPdf {
    PdfEdges edges;
    std::vector<std::array<double, 3>> moments;
    std::vector<std::int64_t> kept;
};

// Slot 0 is the full sample; slot r+1 omits region r.
using ComponentPdfs = std::array<ReplicaPdf, 2>;
using ReplicaCosts = std::vector<std::vector<double>>;
using TrialGrids = std::array<const std::vector<double> *, 2>;

ComponentPdfs prepare_replica_pdf(const PairStore &pairs, const std::vector<int> &labels,
                                  int regions, int bins, MPI_Comm comm);
std::array<ReplicaCosts, 2> replica_chi2(const PairStore &pairs, const std::vector<int> &labels,
                                        const ComponentPdfs &pdf, const TrialGrids &grids, MPI_Comm comm);
std::array<std::vector<double>, 2> fit_replica_pdf(const PairStore &pairs,
                                                const std::vector<int> &labels,
                                                const ComponentPdfs &pdf, MPI_Comm comm);
} // namespace bcg
