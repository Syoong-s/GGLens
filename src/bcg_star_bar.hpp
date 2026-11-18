#pragma once
#include "bcg_stream.hpp"

namespace bcg {
// ==========================================
// Structure: StarBarCuts
// Method: Hold the shared global threshold or compact sorted exposure thresholds.
// ==========================================
struct StarBarCuts {
    bool enabled = false, exposure = false;
    double global = 0;
    std::vector<int> ids;
    std::vector<double> thresholds;
    bool accepts(const BarDatum &b) const;
};
StarBarCuts fit_star_bar(const BarRows &bar, const Cuts &cuts, MPI_Comm comm, std::size_t budget);
void apply_star_bar(LocalSources &local, const Cuts &cuts, MPI_Comm comm);
} // namespace bcg
