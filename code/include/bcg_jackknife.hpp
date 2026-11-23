#pragma once
#include "bcg_pair_store.hpp"

namespace bcg {
using PdfEdges = std::vector<std::vector<double>>;

struct ReplicaPdf {
    PdfEdges edges;
    std::vector<std::array<double, 3>> moments;
    std::vector<std::int64_t> kept;
};

// Slot 0 is the full sample; slot r+1 omits region r.
ReplicaPdf prepare_replica_pdf(const PairStore &pairs, const std::vector<int> &labels,
                                int regions, int component, int bins, MPI_Comm comm);
std::vector<std::vector<double>> replica_chi2(
    const PairStore &pairs, const std::vector<int> &labels, int component,
    const PdfEdges &edges, const std::vector<double> &trials, MPI_Comm comm);
std::vector<double> fit_replica_pdf(const PairStore &pairs, const std::vector<int> &labels,
                                    int component, const ReplicaPdf &pdf, MPI_Comm comm);
} // namespace bcg
