#include "bcg_pair_store.hpp"

namespace bcg {
// ==========================================
// Function: PairStore constructor
// Method: Limit all bin write buffers jointly through the caller's divided quota.
// ==========================================
PairStore::PairStore(const std::string &file, std::size_t read_bytes, std::size_t write_bytes,
                     std::size_t lenses, int regions, bool field_distortion)
    : file_(file), read_bytes_(read_bytes), field_distortion_(field_distortion),
      present(lenses, 0), region_count(regions, 0) {
    if (field_distortion_) fd_writer_ = std::make_unique<BinaryWriter<FdRecord>>(file, write_bytes);
    else writer_ = std::make_unique<BinaryWriter<PairRecord>>(file, write_bytes);
}

// ==========================================
// Function: PairStore::append
// Method: Write each selected pair once and retain only small radial/lens metadata.
// ==========================================
void PairStore::append(const PairRecord &pair, double radius, int region) {
    if (!writer_ && !fd_writer_) throw std::runtime_error("Pair spool is already closed");
    if (count == std::numeric_limits<std::int64_t>::max())
        throw std::runtime_error("Pair count exceeds 64-bit limit");
    if (field_distortion_) fd_writer_->append({pair.gt, pair.det});
    else writer_->append(pair);
    ++count;
    log_radius_sum += std::log(radius);
    present.at(static_cast<std::size_t>(pair.lens)) = 1; ++region_count.at(region);
}

// ==========================================
// Function: PairStore::finish
// Method: Close a bin spool and immediately release its bounded write allocation.
// ==========================================
void PairStore::finish() {
    if (writer_) { writer_->close(); writer_.reset(); }
    if (fd_writer_) { fd_writer_->close(); fd_writer_.reset(); }
}

// ==========================================
// Function: PairStore::load_if_small
// Method: Cache only one bin if its exact record payload fits the resident quota.
// ==========================================
bool PairStore::load_if_small(std::size_t bytes) {
    if (static_cast<std::uint64_t>(count) > bytes / sizeof(PairRecord)) return false;
    if (field_distortion_) {
        memory_.clear();
        memory_.reserve(static_cast<std::size_t>(count));
        if (count) scan_disk([&](const PairRecord &pair) { memory_.push_back(pair); });
        return true;
    }
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
    if (memory_.empty()) scan_disk(visit);
    else for (const auto &p : memory_) visit(p);
}

void PairStore::scan_disk(const std::function<void(const PairRecord &)> &visit) const {
    std::int64_t seen = 0;
    if (field_distortion_) {
        BinaryReader<FdRecord> in(file_, read_bytes_);
        FdRecord record;
        while (in.next(record)) {
            if (seen == count) throw std::runtime_error("Pair spool has extra records: " + file_);
            visit({record[0], 0, record[1], 0, 0, 0});
            ++seen;
        }
    } else {
        BinaryReader<PairRecord> in(file_, read_bytes_);
        PairRecord pair;
        while (in.next(pair)) {
            if (seen == count) throw std::runtime_error("Pair spool has extra records: " + file_);
            visit(pair);
            ++seen;
        }
    }
    if (seen != count) throw std::runtime_error("Pair spool count changed: " + file_);
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
