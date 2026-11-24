#pragma once

#include "bcg_stream.hpp"

namespace bcg {
// ==========================================
// Structure: PairRecord
// Method: Store four double PDF inputs and the lens ID; accumulate R separately.
// ==========================================
struct PairRecord {
    double gt = 0, gx = 0, det = 0, dex = 0;
    std::int32_t lens = 0;
    std::int32_t padding = 0;
};
static_assert(sizeof(PairRecord) == 40, "Compact pair layout must use forty bytes");

// ==========================================
// Class: PairStore
// Method: Spool one R bin; load only a budget-sized bin, otherwise scan disk blocks.
// ==========================================
class PairStore {
    using FdRecord = std::array<double, 2>;
    std::string file_;
    std::size_t read_bytes_;
    bool field_distortion_;
    std::unique_ptr<BinaryWriter<PairRecord>> writer_;
    std::unique_ptr<BinaryWriter<FdRecord>> fd_writer_;
    std::vector<PairRecord> memory_;
    void scan_disk(const std::function<void(const PairRecord &)> &visit) const;
public:
    std::int64_t count = 0;
    double log_radius_sum = 0;
    std::vector<unsigned char> present;
    std::vector<std::int64_t> region_count;
    PairStore(const std::string &file, std::size_t read_bytes, std::size_t write_bytes,
              std::size_t lenses, int regions, bool field_distortion = false);
    void append(const PairRecord &pair, double radius, int region = 0);
    void finish();
    bool load_if_small(std::size_t bytes);
    void scan(const std::function<void(const PairRecord &)> &visit) const;
    void release();
};

} // namespace bcg
