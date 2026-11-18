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
    std::string file_;
    std::size_t read_bytes_;
    std::unique_ptr<BinaryWriter<PairRecord>> writer_;
    std::vector<PairRecord> memory_;
public:
    std::int64_t count = 0;
    double log_radius_sum = 0;
    std::vector<unsigned char> present;
    std::vector<std::int64_t> region_count;
    PairStore(const std::string &file, std::size_t read_bytes, std::size_t write_bytes,
              std::size_t lenses, int regions);
    void append(const PairRecord &pair, double radius, int region);
    void finish();
    bool load_if_small(std::size_t bytes);
    void scan(const std::function<void(const PairRecord &)> &visit) const;
    void release();
};

// ==========================================
// Structure: PairView
// Method: Select a component and omitted region without copying pair measurements.
// ==========================================
struct PairView {
    const PairStore &store;
    const std::vector<int> &labels;
    int component = 0, omitted = -1;
    // ==========================================
    // Function: size
    // Method: Obtain local leave-one-region-out counts from pairing metadata.
    // ==========================================
    std::int64_t size() const {
        return store.count - (omitted < 0 ? 0 : store.region_count.at(omitted));
    }
    // ==========================================
    // Function: scan
    // Method: Feed component values directly to statistics, applying region exclusion.
    // ==========================================
    void scan(const std::function<void(double, double)> &visit) const {
        store.scan([&](const PairRecord &p) {
            if (omitted >= 0 && labels[p.lens] == omitted) return;
            visit(component ? p.gx : p.gt, component ? p.dex : p.det);
        });
    }
};
} // namespace bcg
