#include "bcg_dr9.hpp"
#include "bcg_mpi_rows.hpp"
#include "bcg_dr9_rows.hpp"
#include <iostream>
#include <map>

namespace bcg {
namespace {
struct Member {
    std::size_t group_index;
    std::int64_t galaxy_id;
};
struct Galaxy {
    std::int64_t id;
    double ra, dec, stellar_mass;
    int select;
};
} // namespace

// ==========================================
// Function: select_dr9_groups
// Method: Read group ID, group redshift and selection mass in parallel.
// ==========================================
std::vector<Lens> select_dr9_groups(const CatalogTriple &triple, const GroupColumns &columns,
                                    double zmin, double zmax, double mmin, double mmax, MPI_Comm comm) {
    std::vector<Lens> local;
    const auto rows = read_dr9_range(triple.group, columns.text, comm, [&](Dr9Fields &row, std::int64_t number) {
        Lens lens;
        lens.id = columns.id == 0 ? number : row.number<std::int64_t>(columns.id);
        lens.z = row.number<double>(columns.redshift);
        if (!(lens.z >= zmin && lens.z < zmax)) return;
        const double mass = row.number<double>(columns.halo_mass);
        if (columns.mass_format == MassFormat::linear && !(mass > 0)) return;
        lens.log_mh = columns.mass_format == MassFormat::linear ? std::log10(mass) : mass;
        if (!(lens.log_mh >= mmin && lens.log_mh < mmax)) return;
        if (lens.id <= 0) throw std::runtime_error("Invalid selected Group ID in " + triple.group);
        lens.catalog_id = triple.id;
        local.push_back(lens);
    });
    if (columns.id == 0) {
        const auto offset = preceding_catalog_rows(rows, comm);
        for (auto &lens : local) lens.id += offset;
    }
    return allgather_rows(local, comm);
}

// ==========================================
// Function: match_dr9_members
// Method: Resolve exactly one brightest galaxy per selected group within this block.
// ==========================================
void match_dr9_members(const CatalogTriple &triple, const MemberColumns &columns,
                       std::vector<Lens> &lenses, MPI_Comm comm) {
    if (lenses.empty()) return;
    std::unordered_map<std::int64_t, std::size_t> wanted;
    wanted.reserve(lenses.size());
    for (std::size_t i = 0; i < lenses.size(); ++i) {
        if (!wanted.emplace(lenses[i].id, i).second)
            throw std::runtime_error("Duplicate group ID in " + triple.group);
    }
    std::vector<Member> local;
    read_dr9_range(triple.member, columns.text, comm, [&](Dr9Fields &row, std::int64_t) {
        const auto galaxy_id = row.number<std::int64_t>(columns.galaxy_id);
        const auto group_id = row.number<std::int64_t>(columns.group_id);
        const auto brightest = row.number<int>(columns.brightest);
        if (brightest != columns.brightest_value) return;
        const auto it = wanted.find(group_id);
        if (it != wanted.end()) local.push_back({it->second, galaxy_id});
    });
    for (const auto &member : allgather_rows(local, comm)) {
        auto &lens = lenses[member.group_index];
        if (lens.galaxy_id >= 0 || member.galaxy_id < 0)
            throw std::runtime_error("Invalid or multiple BCGs for group " + std::to_string(lens.id) +
                                     " in " + triple.member);
        lens.galaxy_id = member.galaxy_id;
    }
    for (const auto &lens : lenses)
        if (lens.galaxy_id < 0)
            throw std::runtime_error("Missing BCG for group " + std::to_string(lens.id) +
                                     " in " + triple.member);
}

// ==========================================
// Function: resolve_dr9_galaxies
// Method: Union requested IDs by canonical galaxy file, then scan each unique file once.
// ==========================================
void resolve_dr9_galaxies(const std::vector<CatalogTriple> &triples, const GalaxyColumns &columns,
                          std::vector<Lens> &lenses, MPI_Comm comm) {
    std::unordered_map<int, std::string> path_by_block;
    for (const auto &triple : triples) path_by_block.emplace(triple.id, triple.galaxy);
    std::map<std::string, std::unordered_map<std::int64_t, std::size_t>> wanted_by_path;
    for (std::size_t i = 0; i < lenses.size(); ++i) {
        const auto &lens = lenses[i];
        auto &wanted = wanted_by_path[path_by_block.at(lens.catalog_id)];
        if (!wanted.emplace(lens.galaxy_id, i).second)
            throw std::runtime_error("Repeated BCG ID in shared Galaxy file: " +
                                     std::to_string(lens.galaxy_id));
    }

    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    for (const auto &[path, wanted] : wanted_by_path) {
        if (!rank) std::cout << "Galaxy parallel scan: " << path
                             << " (" << wanted.size() << " unique BCG IDs)" << std::endl;
        const auto offset = galaxy_row_offset(path, columns, comm);
        std::vector<Galaxy> local;
        read_dr9_range(path, columns.text, comm, [&](Dr9Fields &row, std::int64_t number) {
            const auto id = columns.id == 0 ? offset + number : row.number<std::int64_t>(columns.id);
            if (!wanted.count(id)) return;
            Galaxy galaxy{id, row.number<double>(columns.ra), row.number<double>(columns.dec), 0, 0};
            const double mass = row.number<double>(columns.stellar_mass);
            galaxy.stellar_mass = columns.mass_format == MassFormat::log10 ? std::pow(10.0, mass) : mass;
            if (columns.photometry_select != 0) galaxy.select = row.number<int>(columns.photometry_select);
            if (!(galaxy.ra >= 0 && galaxy.ra < 360 && galaxy.dec >= -90 && galaxy.dec <= 90 &&
                  galaxy.stellar_mass > 0 && std::isfinite(galaxy.stellar_mass)))
                throw std::runtime_error("Invalid BCG coordinates or stellar mass for Galaxy ID " +
                                         std::to_string(id) + " in " + path);
            local.push_back(galaxy);
        });
        for (const auto &galaxy : allgather_rows(local, comm)) {
            auto &lens = lenses[wanted.at(galaxy.id)];
            if (lens.stellar_mass > 0)
                throw std::runtime_error("Duplicate Galaxy ID " + std::to_string(galaxy.id) + " in " + path);
            lens.ra = galaxy.ra;
            lens.dec = galaxy.dec;
            lens.stellar_mass = galaxy.stellar_mass;
            lens.photometry_select = galaxy.select;
        }
        for (const auto &[id, index] : wanted)
            if (lenses[index].stellar_mass == 0)
                throw std::runtime_error("Missing Galaxy ID " + std::to_string(id) + " in " + path);
    }
}
} // namespace bcg
