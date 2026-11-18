#include "bcg_esd.hpp"
#include "bcg_physics.hpp"

namespace bcg_esd {
// ==========================================
// Function: initial_regions
// Method: Keep the original BCG spherical k-means labels and MPI broadcast.
// ==========================================
std::vector<int> initial_regions(const std::vector<bcg::Lens> &lenses, int regions, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    std::vector<int> labels;
    if (!rank) labels = bcg::assign_jackknife(lenses, regions);
    else labels.resize(lenses.size());
    MPI_Bcast(labels.data(), static_cast<int>(labels.size()), MPI_INT, 0, comm); return labels;
}
} // namespace bcg_esd
