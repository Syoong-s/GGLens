#include "bcg_dr9_join.hpp"
#include "bcg_dr9_rows.hpp"
#include <iostream>

namespace bcg {
namespace {
struct OwnedRequest {
    std::size_t group;
    bool found = false;
};
struct GalaxyPosition {
    std::int64_t id;
    MemberPosition position;
};

// ==========================================
// Function: owner
// Method: Give an explicit Galaxy ID the same owner during both catalog scans.
// ==========================================
int owner(std::int64_t id, int ranks) {
    return static_cast<int>(id % ranks);
}

// ==========================================
// Class: RequestedIds
// Method: Filter the galaxy scan with an exact compact shared set of selected IDs.
// ==========================================
class RequestedIds {
    std::int64_t first_ = 0, last_ = -1;
    std::vector<std::uint64_t> bitmap_;
    std::unordered_set<std::int64_t> sparse_;
public:
    // ==========================================
    // Function: RequestedIds
    // Method: Use a bounded bitmap for dense spans and compact IDs for small or sparse selections.
    // ==========================================
    RequestedIds(const std::vector<MemberRequest> &local, MPI_Comm comm) {
        first_ = std::numeric_limits<std::int64_t>::max();
        for (const auto &row : local) {
            first_ = std::min(first_, row.galaxy_id);
            last_ = std::max(last_, row.galaxy_id);
        }
        MPI_Allreduce(MPI_IN_PLACE, &first_, 1, MPI_INT64_T, MPI_MIN, comm);
        MPI_Allreduce(MPI_IN_PLACE, &last_, 1, MPI_INT64_T, MPI_MAX, comm);
        if (last_ < first_) return;
        std::uint64_t count = local.size();
        MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_UINT64_T, MPI_SUM, comm);
        const auto words = static_cast<std::uint64_t>(last_ - first_) / 64 + 1;
        constexpr auto max_words = (32U << 20) / sizeof(std::uint64_t);
        if (words <= max_words && words <= 8 * count) {
            bitmap_.resize(static_cast<std::size_t>(words), 0);
            for (const auto &row : local) {
                const auto offset = static_cast<std::uint64_t>(row.galaxy_id - first_);
                bitmap_[offset / 64] |= std::uint64_t{1} << (offset % 64);
            }
            MPI_Allreduce(MPI_IN_PLACE, bitmap_.data(), static_cast<int>(bitmap_.size()),
                          MPI_UINT64_T, MPI_BOR, comm);
        } else {
            std::vector<std::int64_t> ids;
            ids.reserve(local.size());
            for (const auto &row : local) ids.push_back(row.galaxy_id);
            const auto all_ids = allgather_rows(ids, comm);
            sparse_.insert(all_ids.begin(), all_ids.end());
        }
    }

    // ==========================================
    // Function: contains
    // Method: Select galaxy rows by ID regardless of where they appear in the file.
    // ==========================================
    bool contains(std::int64_t id) const {
        if (id < first_ || id > last_) return false;
        if (bitmap_.empty()) return sparse_.count(id) != 0;
        const auto offset = static_cast<std::uint64_t>(id - first_);
        return (bitmap_[offset / 64] >> (offset % 64)) & std::uint64_t{1};
    }

    bool empty() const { return last_ < first_; }
};

// ==========================================
// Function: read_position
// Method: Read coordinates and, for LWC, the z-band luminosity from the selected Galaxy row.
// ==========================================
MemberPosition read_position(Dr9Fields &row, std::int64_t id, bool luminosity_weighted) {
    MemberPosition position{row.number<double>(), row.number<double>(), 1.0};
    if (luminosity_weighted) {
        row.take(); // Redshift.
        row.take(); // z-band absolute magnitude.
        position.weight = row.number<double>();
    }
    if (!(position.ra >= 0 && position.ra < 360 && position.dec >= -90 && position.dec <= 90))
        throw std::runtime_error("Invalid DR9 coordinates for Galaxy ID " + std::to_string(id));
    if (!(position.weight > 0 && std::isfinite(position.weight)))
        throw std::runtime_error("LWC requires positive z-band luminosity for Galaxy ID " + std::to_string(id));
    return position;
}
} // namespace

// ==========================================
// Function: join_dr9_member_positions
// Method: Route requests once, then consume matching rows during a single parallel galaxy scan.
// ==========================================
void join_dr9_member_positions(
    const std::string &path, const std::vector<MemberRequest> &requests, bool luminosity_weighted,
    const std::function<void(std::size_t, const MemberPosition &)> &accumulate, MPI_Comm comm) {
    int rank = 0, ranks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    const RequestedIds ids(requests, comm);
    if (ids.empty()) return;

    std::unordered_map<std::int64_t, OwnedRequest> owned;
    {
        std::vector<std::vector<MemberRequest>> outgoing(ranks);
        for (const auto &req : requests) outgoing[owner(req.galaxy_id, ranks)].push_back(req);
        const auto incoming = exchange_rows(outgoing, comm);
        owned.reserve(incoming.size());
        for (const auto &req : incoming) {
            if (!owned.emplace(req.galaxy_id, OwnedRequest{req.group}).second)
                throw std::runtime_error("Duplicate selected member Galaxy ID " +
                                         std::to_string(req.galaxy_id) + " for " + path);
        }
    }

    if (!rank) std::cout << "DR9 Galaxy parallel scan: " << path << std::endl;
    Dr9Rows input(path, comm);
    std::string line;
    bool exhausted = false;
    std::vector<std::vector<GalaxyPosition>> outgoing(ranks);
    constexpr int batch_rows = 262144;
    for (;;) {
        int selected = 0;
        for (int i = 0; i < batch_rows && !exhausted; ++i) {
            if (!input.next(line)) { exhausted = true; break; }
            Dr9Fields row(line, path);
            const auto id = row.number<std::int64_t>();
            if (!ids.contains(id)) continue;
            const auto position = read_position(row, id, luminosity_weighted);
            outgoing[owner(id, ranks)].push_back({id, position});
            selected = 1;
        }
        int active[] = {!exhausted, selected};
        MPI_Allreduce(MPI_IN_PLACE, active, 2, MPI_INT, MPI_MAX, comm);
        if (active[1]) {
            for (const auto &row : exchange_rows(outgoing, comm)) {
                auto &entry = owned.at(row.id);
                if (entry.found)
                    throw std::runtime_error("Duplicate Galaxy ID " + std::to_string(row.id) + " in " + path);
                accumulate(entry.group, row.position);
                entry.found = true;
            }
            for (auto &bucket : outgoing) bucket.clear();
        }
        if (!active[0]) break;
    }
    for (const auto &[id, entry] : owned)
        if (!entry.found)
            throw std::runtime_error("Missing member Galaxy ID " + std::to_string(id) + " in " + path);
}
} // namespace bcg
