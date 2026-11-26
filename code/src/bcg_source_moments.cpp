#include "bcg_source_moments.hpp"

namespace bcg {
SourceMomentStore::SourceMomentStore(const std::string &file, std::size_t io_bytes, bool self_calibrate,
                                     const Cosmology &cosmology)
    : file_(file), io_bytes_(io_bytes), cosmology_(cosmology) {
    if (self_calibrate) raw_writer_ = std::make_unique<BinaryWriter<RawSource>>(file_, io_bytes_);
    else writer_ = std::make_unique<BinaryWriter<DiskSource>>(file_, io_bytes_);
}

SourceMomentStore::~SourceMomentStore() {
    writer_.reset();
    raw_writer_.reset();
    if (reader_.is_open()) reader_.close();
    std::error_code error;
    fs::remove(file_, error);
}

std::uint32_t SourceMomentStore::block_id() const {
    if (blocks_.size() >= std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("Too many matched source chunks");
    return static_cast<std::uint32_t>(blocks_.size());
}

void SourceMomentStore::append(const Source &source, std::uint32_t order) {
    const DiskSource row{source.g1, source.g2, source.de, source.h1, source.h2,
                         source.ra, source.dec, source.zp, order};
    if (raw_writer_) raw_writer_->append({row, source.gf1, source.gf2});
    else writer_->append(row);
    ++rows_;
}

void SourceMomentStore::finish_block(std::size_t slots) {
    if (rows_ == block_start_) return;
    blocks_.push_back({block_start_, static_cast<std::uint32_t>(rows_ - block_start_),
                       static_cast<std::uint32_t>(slots)});
    block_start_ = rows_;
}

void SourceMomentStore::finish(std::vector<LensGeometry> lenses) {
    if (writer_) {
        writer_->close();
        writer_.reset();
    }
    if (raw_writer_) {
        raw_writer_->close();
        raw_writer_.reset();
    }
    lenses_ = std::move(lenses);
}

void SourceMomentStore::calibrate(const ShearCalibration &calibration) {
    if (!rows_) return;

    // Read ahead before replacing each raw block. The smaller output records
    // always end before the next unread block, so there is no second table.
    const auto block_rows = std::max<std::size_t>(1, io_bytes_ / sizeof(RawSource));
    std::vector<RawSource> raw(std::min<std::uint64_t>(block_rows, rows_));
    std::vector<DiskSource> calibrated(raw.size());
    std::ifstream in(file_, std::ios::binary);
    std::fstream out(file_, std::ios::binary | std::ios::in | std::ios::out);
    for (std::uint64_t remaining = rows_; remaining;) {
        const auto n = static_cast<std::size_t>(std::min<std::uint64_t>(remaining, raw.size()));
        if (!in.read(reinterpret_cast<char *>(raw.data()), n * sizeof(RawSource)))
            throw std::runtime_error("Cannot read raw source moments: " + file_);
        for (std::size_t i = 0; i < n; ++i) {
            const auto &row = raw[i].source;
            Source source;
            source.g1 = row.g1;
            source.g2 = row.g2;
            source.de = row.de;
            source.h1 = row.h1;
            source.h2 = row.h2;
            source.gf1 = raw[i].gf1;
            source.gf2 = raw[i].gf2;
            calibration.apply(source);
            calibrated[i] = row;
            calibrated[i].g1 = source.g1;
            calibrated[i].g2 = source.g2;
        }
        if (!out.write(reinterpret_cast<const char *>(calibrated.data()), n * sizeof(DiskSource)))
            throw std::runtime_error("Cannot write calibrated source moments: " + file_);
        remaining -= n;
    }
    if (in.peek() != std::char_traits<char>::eof() || in.bad())
        throw std::runtime_error("Invalid raw source spool length: " + file_);
    out.close();
    if (!out) throw std::runtime_error("Cannot close calibrated source spool: " + file_);
    in.close();
    fs::resize_file(file_, rows_ * sizeof(DiskSource));
}

const std::vector<PreparedSource> &SourceMomentStore::load_block(std::uint32_t id) {
    if (loaded_ == id) return memory_;
    if (!reader_.is_open()) reader_.open(file_, std::ios::binary);
    const auto &block = blocks_.at(id);
    reader_.seekg(static_cast<std::streamoff>(block.start * sizeof(DiskSource)));
    memory_.resize(block.slots);
    const auto block_rows = std::max<std::size_t>(1, io_bytes_ / sizeof(DiskSource));
    std::vector<DiskSource> rows(std::min<std::size_t>(block_rows, block.count));
    for (std::size_t remaining = block.count; remaining;) {
        const auto n = std::min(remaining, rows.size());
        if (!reader_.read(reinterpret_cast<char *>(rows.data()), n * sizeof(DiskSource)))
            throw std::runtime_error("Cannot load matched source moments: " + file_);
        for (std::size_t i = 0; i < n; ++i) {
            const auto &row = rows[i];
            memory_.at(row.order) = {{row.g1, row.g2, row.de, row.h1, row.h2},
                                      row.zp, cosmology_.chi(row.zp), sky_xyz(row.ra, row.dec)};
        }
        remaining -= n;
    }
    loaded_ = id;
    return memory_;
}
} // namespace bcg
