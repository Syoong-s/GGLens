#pragma once
#include "bcg_catalog.hpp"

namespace bcg {
enum class Center { bcg, lwc, nwc };

// ==========================================
// Function: center_name
// Method: Use the INI spelling in provenance, logs, and cache keys.
// ==========================================
constexpr const char *center_name(Center center) {
    switch (center) {
        case Center::bcg: return "bcg";
        case Center::lwc: return "lwc";
        case Center::nwc: return "nwc";
    }
    throw std::logic_error("Unknown lens center");
}

// ==========================================
// Function: parse_center
// Method: Parse a case-insensitive center selection from the INI.
// ==========================================
inline Center parse_center(std::string name) {
    for (char &letter : name) letter = static_cast<char>(std::tolower(static_cast<unsigned char>(letter)));
    for (const auto center : {Center::bcg, Center::lwc, Center::nwc})
        if (name == center_name(center)) return center;
    throw std::invalid_argument("center must be bcg, lwc or nwc");
}

std::vector<CatalogTriple> read_group_manifest(const std::string &path);
std::string catalog_fingerprint(const std::vector<CatalogTriple> &triples,
                                double zmin, double zmax, double mmin, double mmax, Center center);
bool read_lens_cache(const std::string &path, const std::string &key,
                     const std::vector<CatalogTriple> &triples, Center center,
                     std::vector<Lens> &lenses);
void write_lens_cache(const std::string &path, const std::string &key,
                      const std::vector<Lens> &lenses);

std::vector<Lens> select_dr9_groups(const CatalogTriple &triple, double zmin, double zmax,
                                    double mmin, double mmax, MPI_Comm comm);
void match_dr9_members(const CatalogTriple &triple, std::vector<Lens> &lenses, MPI_Comm comm);
void resolve_dr9_galaxies(const std::vector<CatalogTriple> &triples,
                          std::vector<Lens> &lenses, MPI_Comm comm);
void resolve_dr9_centroids(const std::vector<CatalogTriple> &triples,
                           std::vector<Lens> &lenses, Center center, MPI_Comm comm);
} // namespace bcg
