#include "bcg_dr9.hpp"
#include "bcg_esd.hpp"
#include "bcg_mpi_rows.hpp"
#include <iostream>

namespace bcg_esd {
// ==========================================
// Function: load_lenses
// Method: Reuse a compatible cache or match the numbered DR9 catalogs on all ranks.
// ==========================================
std::vector<bcg::Lens> load_lenses(const Config &cfg, MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    std::vector<bcg::Lens> lenses;
    std::vector<bcg::CatalogTriple> triples;
    std::vector<std::string> manifest;
    std::string key;
    int cache_hit = 0;
    if (!rank) {
        triples = bcg::read_group_manifest(cfg.group_list);
        if (!cfg.cache.empty()) {
            key = bcg::catalog_fingerprint(triples, cfg.zmin, cfg.zmax, cfg.mmin, cfg.mmax, cfg.center);
            cache_hit = bcg::read_lens_cache(cfg.cache, key, triples, cfg.center, lenses);
        }
        for (const auto &triple : triples) {
            manifest.insert(manifest.end(), {std::to_string(triple.id), triple.group,
                                             triple.galaxy, triple.member});
        }
    }
    MPI_Bcast(&cache_hit, 1, MPI_INT, 0, comm);
    if (cache_hit) {
        lenses = bcg::allgather_rows(lenses, comm);
        if (!rank) std::cout << "Loaded " << lenses.size() << " DR9 " << bcg::center_name(cfg.center)
                             << " centers from cache " << cfg.cache << std::endl;
        return lenses;
    }
    bcg::broadcast_paths(manifest, comm);
    if (rank) {
        for (std::size_t i = 0; i < manifest.size(); i += 4)
            triples.push_back({std::stoi(manifest[i]), manifest[i + 1], manifest[i + 2], manifest[i + 3]});
    }
    for (const auto &triple : triples) {
        auto selected = bcg::select_dr9_groups(triple, cfg.zmin, cfg.zmax, cfg.mmin, cfg.mmax, comm);
        if (cfg.center == bcg::Center::bcg) bcg::match_dr9_members(triple, selected, comm);
        if (!rank) std::cout << "DR9 block " << triple.id << ": " << selected.size() << " selected groups\n";
        lenses.insert(lenses.end(), selected.begin(), selected.end());
    }
    if (lenses.empty()) throw std::runtime_error("No DR9 groups pass the redshift and mass cuts");
    if (cfg.center == bcg::Center::bcg) bcg::resolve_dr9_galaxies(triples, lenses, comm);
    else bcg::resolve_dr9_centroids(triples, lenses, cfg.center, comm);
    if (!rank) {
        bcg::write_lens_cache(cfg.cache, key, lenses);
        std::cout << "Resolved " << lenses.size() << " DR9 " << bcg::center_name(cfg.center)
                  << " centers" << std::endl;
    }
    return lenses;
}
} // namespace bcg_esd
