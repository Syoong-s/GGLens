#include "bcg_pair_store.hpp"

namespace bcg {
namespace {
template<class Record>
void load_records(const std::string &file, std::vector<Record> &records, std::size_t count) {
    records.resize(count);
    if (records.empty()) return;
    std::ifstream in(file, std::ios::binary);
    in.read(reinterpret_cast<char *>(records.data()), records.size() * sizeof(Record));
    if (!in || in.peek() != std::char_traits<char>::eof() || in.bad())
        throw std::runtime_error("Invalid pair spool length: " + file);
}

// A nonzero varint advances the source index. Zero introduces a tag.
// An even tag selects a lens. An odd tag gives a relative base for a 64-bit bitmap.
void write_varint(BinaryWriter<unsigned char> &out, std::uint32_t value) {
    while (value >= 128) {
        out.append(static_cast<unsigned char>(128 | (value & 127)));
        value >>= 7;
    }
    out.append(static_cast<unsigned char>(value));
}

std::uint32_t read_varint(BinaryReader<unsigned char> &in) {
    unsigned char byte;
    if (!in.next(byte)) throw std::runtime_error("Truncated compressed pair stream");
    if (byte < 128) return byte;
    std::uint32_t result = byte & 127;
    for (int shift = 7; shift <= 28; shift += 7) {
        if (!in.next(byte)) throw std::runtime_error("Truncated compressed pair stream");
        if (shift == 28 && byte > 15) throw std::runtime_error("Overlong compressed pair index");
        result |= static_cast<std::uint32_t>(byte & 127) << shift;
        if (byte < 128) return result;
    }
    throw std::runtime_error("Overlong compressed pair index");
}

PairRecord project_geometry(const PreparedSource &s, const LensGeometry &l, std::int32_t lens) {
    const double xy = s.xyz[0] * l.xyz[0] + s.xyz[1] * l.xyz[1];
    const double east = s.xyz[0] * l.xyz[1] - s.xyz[1] * l.xyz[0];
    const double phi = source_bearing(s.xyz.data(), l.xyz, xy, east);
    const double beta = Cosmology::beta(l.z, s.zp, l.chi, s.chi);
    return project_pair(s.moments, phi, beta, lens);
}
} // namespace

// ==========================================
// Function: PairStore constructor
// Method: Limit all bin write buffers jointly through the caller's divided quota.
// ==========================================
PairStore::PairStore(const std::string &file, std::size_t read_bytes, std::size_t write_bytes,
                     std::size_t lenses, int regions, bool field_distortion,
                     std::shared_ptr<SourceMomentStore> sources)
    : file_(file), read_bytes_(read_bytes), field_distortion_(field_distortion), sources_(std::move(sources)),
      present(lens_mask_words(lenses), 0), region_count(regions, 0) {
    if (field_distortion_) fd_writer_ = std::make_unique<BinaryWriter<FdRecord>>(file, write_bytes);
    else if (sources_) link_writer_ = std::make_unique<BinaryWriter<unsigned char>>(file, write_bytes);
    else writer_ = std::make_unique<BinaryWriter<PairRecord>>(file, write_bytes);
}

// ==========================================
// Function: PairStore::append
// Method: Write each selected pair once and retain only small radial/lens metadata.
// ==========================================
void PairStore::append(const PairRecord &pair, double radius, int region) {
    if (!writer_ && !fd_writer_) throw std::runtime_error("Pair spool is already closed");
    note_pair(pair.lens, radius, region);
    if (field_distortion_) fd_writer_->append({pair.gt, pair.det});
    else writer_->append(pair);
}

void PairStore::flush_members() {
    if (!pending_mask_) return;
    // A dense word amortizes the eight bitmap bytes and its tag across the hits.
    if (__builtin_popcountll(pending_mask_) >= 12) {
        write_varint(*link_writer_, 0);
        write_varint(*link_writer_, ((pending_base_ - previous_source_) << 1) | 1);
        for (int shift = 0; shift < 64; shift += 8)
            link_writer_->append(static_cast<unsigned char>(pending_mask_ >> shift));
        previous_source_ = pending_base_ + 64 - __builtin_clzll(pending_mask_);
    } else {
        while (pending_mask_) {
            const auto source = pending_base_ + __builtin_ctzll(pending_mask_) + 1;
            write_varint(*link_writer_, source - previous_source_);
            previous_source_ = source;
            pending_mask_ &= pending_mask_ - 1;
        }
    }
    pending_mask_ = 0;
}

void PairStore::append_link(std::uint32_t block, std::uint32_t source, std::int32_t lens, double radius) {
    if (!link_writer_) throw std::runtime_error("Pair link spool is already closed");
    const bool new_block = link_blocks_.empty() || link_blocks_.back().source_block != block;
    if (new_block || lens != previous_lens_) {
        flush_members();
        if (new_block) link_blocks_.push_back({count, block});
        if (lens < 0) throw std::runtime_error("Invalid indexed lens");
        write_varint(*link_writer_, 0);
        write_varint(*link_writer_, static_cast<std::uint32_t>(lens) << 1);
        previous_lens_ = lens;
        previous_source_ = 0;
    }
    const auto base = source & ~std::uint32_t{63};
    if (base != pending_base_) {
        flush_members();
        pending_base_ = base;
    }
    const auto bit = std::uint64_t{1} << (source & 63);
    if (source < previous_source_ || bit <= pending_mask_)
        throw std::runtime_error("Pair indices must be in KD order");
    pending_mask_ |= bit;
    note_pair(lens, radius, 0);
}

void PairStore::note_pair(std::int32_t lens, double radius, int region) {
    if (count == std::numeric_limits<std::int64_t>::max())
        throw std::runtime_error("Pair count exceeds 64-bit limit");
    ++count;
    log_radius_sum += std::log(radius);
    const auto i = static_cast<std::size_t>(lens);
    present.at(i / lenses_per_word) |= std::uint64_t{1} << (i % lenses_per_word);
    ++region_count.at(region);
}

// ==========================================
// Function: PairStore::finish
// Method: Close a bin spool and immediately release its bounded write allocation.
// ==========================================
void PairStore::finish() {
    if (writer_) { writer_->close(); writer_.reset(); }
    if (fd_writer_) { fd_writer_->close(); fd_writer_.reset(); }
    if (link_writer_) {
        flush_members();
        link_writer_->close();
        link_writer_.reset();
    }
    if (!count) sources_.reset();
}

// ==========================================
// Function: PairStore::prepare_scan
// Method: Cache all or part of an indexed bin during the first scan.
// ==========================================
void PairStore::prepare_scan(std::size_t cache_bytes) {
    const auto record_bytes = field_distortion_ ? sizeof(FdRecord) : sizeof(PairRecord);
    const auto limit = cache_bytes / record_bytes;
    if (sources_) {
        cache_limit_ = std::min<std::uint64_t>(count, limit);
        memory_.reserve(cache_limit_);
    } else if (static_cast<std::uint64_t>(count) <= limit && memory_.empty() && fd_memory_.empty()) {
        if (field_distortion_) load_records(file_, fd_memory_, static_cast<std::size_t>(count));
        else load_records(file_, memory_, static_cast<std::size_t>(count));
    }
}

void PairStore::scan_links(const std::function<void(const PairRecord &)> &visit) const {
    BinaryReader<unsigned char> in(file_, read_bytes_);
    std::uint64_t seen = 0;
    for (std::size_t block = 0; block < link_blocks_.size(); ++block) {
        const auto &sources = sources_->load_block(link_blocks_[block].source_block);
        const auto end = static_cast<std::uint64_t>(
            block + 1 == link_blocks_.size() ? count : link_blocks_[block + 1].first);
        std::uint32_t previous = 0, bitmap_base = 0;
        std::uint64_t bitmap = 0;
        std::int32_t lens = -1;
        const LensGeometry *geometry = nullptr;
        for (; seen < end; ++seen) {
            if (!bitmap) {
                auto delta = read_varint(in);
                if (!delta) {
                    auto tag = read_varint(in);
                    if (!(tag & 1)) {
                        lens = static_cast<std::int32_t>(tag >> 1);
                        geometry = &sources_->lens(static_cast<std::uint32_t>(lens));
                        previous = 0;
                        delta = read_varint(in);
                        if (!delta) tag = read_varint(in);
                    }
                    if (!delta) {
                        if (!(tag & 1)) throw std::runtime_error("Invalid source bitmask tag");
                        if ((tag >> 1) >= sources.size() - previous)
                            throw std::runtime_error("Invalid source bitmap base");
                        bitmap_base = previous + (tag >> 1);
                        for (int shift = 0; shift < 64; shift += 8) {
                            unsigned char byte;
                            if (!in.next(byte)) throw std::runtime_error("Truncated source bitmask");
                            bitmap |= static_cast<std::uint64_t>(byte) << shift;
                        }
                        if (!bitmap) throw std::runtime_error("Empty source bitmask");
                    }
                }
                if (delta > sources.size() - previous)
                    throw std::runtime_error("Invalid source index delta");
                previous += delta;
            }
            if (bitmap) {
                previous = bitmap_base + __builtin_ctzll(bitmap) + 1;
                bitmap &= bitmap - 1;
            }
            if (!geometry || previous == 0 || previous > sources.size())
                throw std::runtime_error("Invalid source index in compact pair stream");
            if (seen < memory_.size()) {
                visit(memory_[seen]);
            } else {
                const auto projected = project_geometry(sources[previous - 1], *geometry, lens);
                if (seen < cache_limit_) memory_.push_back(projected);
                visit(projected);
            }
        }
        if (bitmap) throw std::runtime_error("Source bitmask exceeds chunk pair count");
    }
    unsigned char tail;
    if (in.next(tail)) throw std::runtime_error("Extra data in compressed pair stream");
}

// ==========================================
// Function: PairStore::scan
// Method: Scan one resident bin or reuse a fixed-size disk block without cloning it.
// ==========================================
void PairStore::scan(const std::function<void(const PairRecord &)> &visit) const {
    if (!count) return;
    if (sources_ && memory_.size() != static_cast<std::uint64_t>(count)) {
        scan_links(visit);
        if (memory_.size() == static_cast<std::uint64_t>(count)) fs::remove(file_);
    } else if (!memory_.empty()) for (const auto &p : memory_) visit(p);
    else if (!fd_memory_.empty()) for (const auto &p : fd_memory_) visit({p[0], 0, p[1], 0, 0});
    else scan_disk(visit);
}

void PairStore::scan_disk(const std::function<void(const PairRecord &)> &visit) const {
    std::int64_t seen = 0;
    if (field_distortion_) {
        BinaryReader<FdRecord> in(file_, read_bytes_);
        FdRecord record;
        while (in.next(record)) {
            if (seen == count) throw std::runtime_error("Pair spool has extra records: " + file_);
            visit({record[0], 0, record[1], 0, 0});
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
    finish();
    std::vector<PairRecord>().swap(memory_);
    std::vector<FdRecord>().swap(fd_memory_);
    std::vector<LinkBlock>().swap(link_blocks_);
    std::vector<std::uint64_t>().swap(present);
    std::vector<std::int64_t>().swap(region_count);
    sources_.reset();
    fs::remove(file_);
}
} // namespace bcg
