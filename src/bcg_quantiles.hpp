#pragma once
#include "bcg_pair_store.hpp"

namespace bcg {
// ==========================================
// Structure: InitialQuantiles
// Method: Exact rank-local starting boundaries for full and omitted-region samples.
// ==========================================
struct InitialQuantiles {
    std::vector<std::vector<double>> values; // slot 0 full; slot r+1 excludes region r.
};
InitialQuantiles prepare_quantiles(const PairStore &store, const std::vector<int> &labels,
                                  int component, int boundaries, int regions,
                                  const ScratchSpace &scratch, std::size_t sort_bytes);
} // namespace bcg
