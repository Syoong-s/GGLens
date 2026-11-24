#pragma once
#include "bcg_catalog.hpp"

namespace bcg {
enum class Center { bcg, lwc, nwc };
enum class CatalogDelimiter : char { whitespace = '\0', comma = ',', tab = '\t' };
enum class MassFormat { log10, linear };

struct CatalogTextFormat {
    CatalogDelimiter delimiter = CatalogDelimiter::whitespace;
    bool header = false;
};
struct GroupColumns {
    int id, redshift, halo_mass;
    CatalogTextFormat text;
    MassFormat mass_format = MassFormat::log10;
};
struct GalaxyColumns {
    int id, ra, dec, luminosity, stellar_mass, photometry_select;
    CatalogTextFormat text;
    MassFormat mass_format = MassFormat::log10;
};
struct MemberColumns {
    int galaxy_id, group_id, brightest;
    CatalogTextFormat text;
    int brightest_value = 1;
};
struct LensColumns {
    GroupColumns group;
    GalaxyColumns galaxy;
    MemberColumns member;
};

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
LensColumns read_lens_columns(const std::vector<CatalogTriple> &triples,
                              LensColumns columns, Center center, MPI_Comm comm);
std::string catalog_fingerprint(const std::vector<CatalogTriple> &triples,
                                double zmin, double zmax, double mmin, double mmax, Center center,
                                const LensColumns &columns);
bool read_lens_cache(const std::string &path, const std::string &key,
                     const std::vector<CatalogTriple> &triples, Center center,
                     std::vector<Lens> &lenses);
void write_lens_cache(const std::string &path, const std::string &key, Center center,
                      const std::vector<Lens> &lenses);

std::vector<Lens> select_dr9_groups(const CatalogTriple &triple, const GroupColumns &columns,
                                    double zmin, double zmax, double mmin, double mmax, MPI_Comm comm);
void match_dr9_members(const CatalogTriple &triple, const MemberColumns &columns,
                       std::vector<Lens> &lenses, MPI_Comm comm);
void resolve_dr9_galaxies(const std::vector<CatalogTriple> &triples, const GalaxyColumns &columns,
                          std::vector<Lens> &lenses, MPI_Comm comm);
void resolve_dr9_centroids(const std::vector<CatalogTriple> &triples, const LensColumns &columns,
                           std::vector<Lens> &lenses, MPI_Comm comm);
} // namespace bcg
