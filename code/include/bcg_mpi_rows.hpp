#pragma once
#include "bcg_catalog.hpp"

namespace bcg {
// ==========================================
// Function: partition
// Method: Return a contiguous range with at most one item of size imbalance.
// ==========================================
inline std::pair<std::uint64_t, std::uint64_t> partition(std::uint64_t n, int rank, int ranks) {
    const auto base = n / ranks, extra = n % ranks;
    const auto begin = base * rank + std::min<std::uint64_t>(rank, extra);
    return {begin, begin + base + (static_cast<std::uint64_t>(rank) < extra)};
}

// ==========================================
// Function: allgather_rows
// Method: Exchange only compact selected records, preserving rank and file order.
// ==========================================
template<class Row>
std::vector<Row> allgather_rows(const std::vector<Row> &local, MPI_Comm comm) {
    static_assert(std::is_trivially_copyable_v<Row>);
    int ranks = 0;
    MPI_Comm_size(comm, &ranks);
    const std::uint64_t local_bytes = local.size() * sizeof(Row);
    std::vector<std::uint64_t> bytes(ranks);
    MPI_Allgather(&local_bytes, 1, MPI_UINT64_T, bytes.data(), 1, MPI_UINT64_T, comm);
    std::vector<int> counts(ranks), offsets(ranks);
    std::uint64_t total = 0;
    for (int r = 0; r < ranks; ++r) {
        if (bytes[r] > std::numeric_limits<int>::max() - total)
            throw std::runtime_error("Selected DR9 records exceed MPI message limit");
        offsets[r] = static_cast<int>(total);
        counts[r] = static_cast<int>(bytes[r]);
        total += bytes[r];
    }
    std::vector<Row> rows(total / sizeof(Row));
    MPI_Allgatherv(local.data(), static_cast<int>(local_bytes), MPI_BYTE, rows.data(),
                   counts.data(), offsets.data(), MPI_BYTE, comm);
    return rows;
}
// ==========================================
// Function: exchange_rows
// Method: Route each bucket to its destination rank using counts of complete records.
// ==========================================
template<class Row>
std::vector<Row> exchange_rows(const std::vector<std::vector<Row>> &buckets, MPI_Comm comm) {
    static_assert(std::is_trivially_copyable_v<Row>);
    int ranks = 1;
    MPI_Comm_size(comm, &ranks);
    std::vector<int> send_count(ranks), receive_count(ranks), send_offset(ranks), receive_offset(ranks);
    std::size_t total_send = 0, total_receive = 0;
    const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
    for (int dest = 0; dest < ranks; ++dest) {
        if (buckets[dest].size() > limit - total_send)
            throw std::runtime_error("DR9 member exchange exceeds the MPI send count limit");
        send_offset[dest] = static_cast<int>(total_send);
        send_count[dest] = static_cast<int>(buckets[dest].size());
        total_send += buckets[dest].size();
    }
    MPI_Alltoall(send_count.data(), 1, MPI_INT, receive_count.data(), 1, MPI_INT, comm);
    for (int source = 0; source < ranks; ++source) {
        const auto count = static_cast<std::size_t>(receive_count[source]);
        if (count > limit - total_receive)
            throw std::runtime_error("DR9 member exchange exceeds the MPI receive count limit");
        receive_offset[source] = static_cast<int>(total_receive);
        total_receive += count;
    }

    std::vector<Row> outgoing, incoming(total_receive);
    outgoing.reserve(total_send);
    for (const auto &bucket : buckets) outgoing.insert(outgoing.end(), bucket.begin(), bucket.end());
    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Row), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(outgoing.data(), send_count.data(), send_offset.data(), type,
                  incoming.data(), receive_count.data(), receive_offset.data(), type, comm);
    MPI_Type_free(&type);
    return incoming;
}
} // namespace bcg
