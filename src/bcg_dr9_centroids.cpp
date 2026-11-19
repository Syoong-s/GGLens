#include "bcg_dr9.hpp"
#include "bcg_dr9_join.hpp"
#include "bcg_dr9_rows.hpp"
#include "bcg_physics.hpp"
#include <map>

namespace bcg {
namespace {
// ==========================================
// Function: read_member_requests
// Method: Resolve Group IDs within one catalog block and retain all selected memberships.
// ==========================================
void read_member_requests(const CatalogTriple &triple, const std::vector<Lens> &lenses,
                          std::vector<MemberRequest> &requests, MPI_Comm comm) {
    std::unordered_map<std::int64_t, std::size_t> groups;
    for (std::size_t i = 0; i < lenses.size(); ++i) {
        if (lenses[i].catalog_id != triple.id) continue;
        if (!groups.emplace(lenses[i].id, i).second)
            throw std::runtime_error("Duplicate group ID in " + triple.group);
    }
    if (groups.empty()) return;

    read_dr9_range(triple.member, comm, [&](Dr9Fields row) {
        const auto galaxy_id = row.number<std::int64_t>();
        const auto group_id = row.number<std::int64_t>();
        const auto it = groups.find(group_id);
        if (it == groups.end()) return;
        if (galaxy_id < 0)
            throw std::runtime_error("Invalid Galaxy ID for group " + std::to_string(group_id) +
                                     " in " + triple.member);
        requests.push_back({galaxy_id, it->second});
    });
}

// ==========================================
// Function: assign_centroid
// Method: Convert a summed sky vector to RA and Dec, including the RA wrap.
// ==========================================
void assign_centroid(Lens &lens, const std::array<double, 4> &sum) {
    const auto [x, y, z, weight] = sum;
    const double norm = std::hypot(x, y, z);
    if (!(weight > 0 && std::isfinite(weight) && norm > 1e-12 * weight && std::isfinite(norm)))
        throw std::runtime_error("Undefined DR9 centroid for group " + std::to_string(lens.id) +
                                 " in block " + std::to_string(lens.catalog_id));
    lens.ra = std::atan2(y, x) / deg;
    if (lens.ra < 0) lens.ra += 360;
    if (lens.ra >= 360) lens.ra = 0;
    lens.dec = std::atan2(z, std::hypot(x, y)) / deg;
}
} // namespace

// ==========================================
// Function: resolve_dr9_centroids
// Method: Reduce spherical member sums across ranks and resolve one center per selected group.
// ==========================================
void resolve_dr9_centroids(const std::vector<CatalogTriple> &triples,
                           std::vector<Lens> &lenses, Center center, MPI_Comm comm) {
    std::map<std::string, std::vector<const CatalogTriple *>> by_galaxy_file;
    for (const auto &triple : triples) by_galaxy_file[triple.galaxy].push_back(&triple);

    std::vector<std::array<double, 4>> sums(lenses.size());
    if (sums.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) / 4)
        throw std::runtime_error("DR9 centroid sums exceed MPI count limits");
    for (const auto &[path, blocks] : by_galaxy_file) {
        std::vector<MemberRequest> requests;
        for (const auto *triple : blocks) read_member_requests(*triple, lenses, requests, comm);
        join_dr9_member_positions(path, requests, center == Center::lwc,
            [&](std::size_t i, const MemberPosition &member) {
                const auto position = sky_xyz(member.ra, member.dec);
                for (int j = 0; j < 3; ++j) sums[i][j] += member.weight * position[j];
                sums[i][3] += member.weight;
            }, comm);
    }
    MPI_Allreduce(MPI_IN_PLACE, sums.data(), static_cast<int>(4 * sums.size()), MPI_DOUBLE, MPI_SUM, comm);
    for (std::size_t i = 0; i < lenses.size(); ++i) assign_centroid(lenses[i], sums[i]);
}
} // namespace bcg
