#include "bcg_pair_store.hpp"

namespace bcg {
// ==========================================
// Function: PairStore constructor
// Method: Limit all bin write buffers jointly through the caller's divided quota.
// ==========================================
PairStore::PairStore(const std::string &file, std::size_t read_bytes, std::size_t write_bytes,
                     std::size_t lenses, int regions)
    : file_(file), read_bytes_(read_bytes),
      writer_(new BinaryWriter<PairRecord>(file, write_bytes)),
      present(lenses, 0), region_count(regions, 0) {}

// ==========================================
// Function: PairStore::append
// Method: Write each selected pair once and retain only small radial/lens metadata.
// ==========================================
void PairStore::append(const PairRecord &pair, double radius, int region) {
    if (!writer_) throw std::runtime_error("Pair spool is already closed");
    if (count == std::numeric_limits<std::int64_t>::max())
        throw std::runtime_error("Pair count exceeds 64-bit limit");
    writer_->append(pair); ++count; log_radius_sum += std::log(radius);
    present.at(static_cast<std::size_t>(pair.lens)) = 1; ++region_count.at(region);
}

// ==========================================
// Function: PairStore::finish
// Method: Close a bin spool and immediately release its bounded write allocation.
// ==========================================
void PairStore::finish() {
    if (writer_) { writer_->close(); writer_.reset(); }
}

// ==========================================
// Function: PairStore::load_if_small
// Method: Cache only one bin if its exact record payload fits the resident quota.
// ==========================================
bool PairStore::load_if_small(std::size_t bytes) {
    if (static_cast<std::uint64_t>(count) > bytes / sizeof(PairRecord)) return false;
    memory_.resize(static_cast<std::size_t>(count));
    if (!count) return true;
    std::ifstream in(file_, std::ios::binary);
    in.read(reinterpret_cast<char *>(memory_.data()),
            static_cast<std::streamsize>(memory_.size() * sizeof(PairRecord)));
    if (!in || in.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("Invalid pair spool length: " + file_);
    return true;
}

// ==========================================
// Function: PairStore::scan
// Method: Scan one resident bin or reuse a fixed-size disk block without cloning it.
// ==========================================
void PairStore::scan(const std::function<void(const PairRecord &)> &visit) const {
    if (!count) return;
    if (!memory_.empty()) {
        for (const auto &p : memory_) visit(p);
    } else {
        BinaryReader<PairRecord> in(file_, read_bytes_);
        PairRecord p; std::int64_t seen = 0;
        while (in.next(p)) { visit(p); ++seen; }
        if (seen != count) throw std::runtime_error("Pair spool count changed: " + file_);
    }
}

// ==========================================
// Function: PairStore::release
// Method: Delete each completed radial bin and release metadata and resident records.
// ==========================================
void PairStore::release() {
    finish(); std::vector<PairRecord>().swap(memory_);
    std::vector<unsigned char>().swap(present);
    std::vector<std::int64_t>().swap(region_count);
    fs::remove(file_);
}
} // namespace bcg
