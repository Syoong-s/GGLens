#pragma once

#include "bcg_source_moments.hpp"

namespace bcg {
// ==========================================
// Structure: PairRecord
// Method: Store four double PDF inputs and the lens ID; accumulate R separately.
// ==========================================
#pragma pack(push, 4)
struct PairRecord {
    double gt = 0, gx = 0, det = 0, dex = 0;
    std::int32_t lens = 0;
};
#pragma pack(pop)
static_assert(sizeof(PairRecord) == 36, "A projected pair stores four doubles and a lens index");

inline PairRecord project_pair(const SourceMoments &s, double phi, double beta, std::int32_t lens) {
    const double c2 = std::cos(2 * phi), s2 = std::sin(2 * phi);
    const double c4 = std::cos(4 * phi), s4 = std::sin(4 * phi);
    const double gt = -(s.g1 * c2 + s.g2 * s2);
    const double gx = s.g1 * s2 - s.g2 * c2;
    const double ut = -(s.h1 * c4 + s.h2 * s4);
    return {gt, gx, beta * (s.de + ut), beta * (s.de - ut), lens};
}

constexpr std::size_t lenses_per_word = std::numeric_limits<std::uint64_t>::digits;
inline std::size_t lens_mask_words(std::size_t lenses) {
    return (lenses + lenses_per_word - 1) / lenses_per_word;
}

// ==========================================
// Class: PairStore
// Method: Spool one R bin and stream its indexed pairs with a bounded projection cache.
// ==========================================
class PairStore {
    using FdRecord = std::array<double, 2>;
    struct LinkBlock {
        std::int64_t first;
        std::uint32_t source_block;
    };
    std::string file_;
    std::size_t read_bytes_;
    bool field_distortion_;
    std::unique_ptr<BinaryWriter<PairRecord>> writer_;
    std::unique_ptr<BinaryWriter<FdRecord>> fd_writer_;
    std::unique_ptr<BinaryWriter<unsigned char>> link_writer_;
    std::shared_ptr<SourceMomentStore> sources_;
    std::vector<LinkBlock> link_blocks_;
    std::int32_t previous_lens_ = -1;
    std::uint32_t previous_source_ = 0, pending_base_ = 0;
    std::uint64_t pending_mask_ = 0;
    std::size_t cache_limit_ = 0;
    mutable std::vector<PairRecord> memory_;
    std::vector<FdRecord> fd_memory_;
    void note_pair(std::int32_t lens, double radius, int region);
    void flush_members();
    void scan_links(const std::function<void(const PairRecord &)> &visit) const;
    void scan_disk(const std::function<void(const PairRecord &)> &visit) const;
public:
    std::int64_t count = 0;
    double log_radius_sum = 0;
    std::vector<std::uint64_t> present;
    std::vector<std::int64_t> region_count;
    PairStore(const std::string &file, std::size_t read_bytes, std::size_t write_bytes,
              std::size_t lenses, int regions, bool field_distortion = false,
              std::shared_ptr<SourceMomentStore> sources = {});
    void append(const PairRecord &pair, double radius, int region = 0);
    void append_link(std::uint32_t block, std::uint32_t source, std::int32_t lens, double radius);
    void finish();
    bool indexed() const { return static_cast<bool>(sources_); }
    void prepare_scan(std::size_t cache_bytes);
    void scan(const std::function<void(const PairRecord &)> &visit) const;
    void release();
};

} // namespace bcg
