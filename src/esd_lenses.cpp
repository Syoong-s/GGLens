#include "bcg_esd.hpp"
#include <iostream>

namespace bcg_esd {
// ==========================================
// Function: load_lenses
// Method: Preserve rank-zero list/cache/BCG resolution and broadcast the selected table.
// ==========================================
std::vector<bcg::Lens> load_lenses(const Config &cfg, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    std::vector<bcg::Lens> lenses;
    if (!rank)
    {
        const auto triples = bcg::read_group_manifest(cfg.group_list);
        const auto key = bcg::catalog_fingerprint(triples, cfg.zmin, cfg.zmax, cfg.mmin, cfg.mmax);
        if (bcg::read_lens_cache(cfg.cache, key, triples, lenses))
            std::cout << "Loaded " << lenses.size() << " BCG lenses from cache\n";
        else {
            lenses = bcg::select_group_lenses(triples, cfg.zmin, cfg.zmax, cfg.mmin, cfg.mmax);
            bcg::assign_bcg_coordinates(lenses, triples); bcg::write_lens_cache(cfg.cache, key, lenses);
            std::cout << "Selected " << lenses.size() << " groups and resolved BCG centers\n";
        }
    }
    bcg::broadcast_lenses(lenses, comm); return lenses;
}
} // namespace bcg_esd
