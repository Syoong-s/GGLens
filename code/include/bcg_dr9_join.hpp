#pragma once
#include "bcg_dr9.hpp"

namespace bcg {
struct MemberRequest {
    std::int64_t galaxy_id;
    std::size_t group;
};

struct MemberPosition {
    double ra, dec, weight;
};

// Send matched member positions to their request owners in bounded batches.
void join_dr9_member_positions(
    const std::string &path, const std::vector<MemberRequest> &requests, const GalaxyColumns &columns,
    const std::function<void(std::size_t, const MemberPosition &)> &accumulate, MPI_Comm comm);
} // namespace bcg
