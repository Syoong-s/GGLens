#pragma once

#include "bcg_catalog.hpp"
#include <chrono>
#include <functional>
#include <memory>

namespace bcg {
// ==========================================
// Class: ScratchSpace
// Method: Own a unique MPI-run/rank directory and remove only owned files.
// ==========================================
class ScratchSpace {
    fs::path root_, rank_dir_;
    bool owned_ = false;
public:
    ScratchSpace(const std::string &base, MPI_Comm comm);
    ~ScratchSpace();
    ScratchSpace(const ScratchSpace &) = delete;
    ScratchSpace &operator=(const ScratchSpace &) = delete;
    std::string path(const std::string &name) const;
    void cleanup(MPI_Comm comm);
};

// ==========================================
// Class: BinaryWriter
// Method: Write POD records through a fixed-capacity buffer, checking every flush.
// ==========================================
template<class T> class BinaryWriter {
    static_assert(std::is_trivially_copyable<T>::value, "Binary records must be POD");
    std::string path_;
    std::ofstream out_;
    std::vector<T> buffer_;
    std::size_t capacity_;
public:
    // ==========================================
    // Function: BinaryWriter constructor
    // Method: Reserve a bounded record buffer before opening a scratch file.
    // ==========================================
    BinaryWriter(const std::string &path, std::size_t bytes)
        : path_(path), capacity_(std::max<std::size_t>(bytes / sizeof(T), 1)) {
        buffer_.reserve(capacity_);
        out_.open(path, std::ios::binary | std::ios::trunc);
        if (!out_) throw std::runtime_error("Cannot create scratch file: " + path);
    }
    // ==========================================
    // Function: append
    // Method: Flush before the fixed-capacity buffer could grow.
    // ==========================================
    void append(const T &value) {
        if (buffer_.size() == capacity_) flush();
        buffer_.push_back(value);
    }
    // ==========================================
    // Function: flush
    // Method: Detect disk-full and incomplete binary writes immediately.
    // ==========================================
    void flush() {
        if (buffer_.empty()) return;
        out_.write(reinterpret_cast<const char *>(buffer_.data()),
                   static_cast<std::streamsize>(buffer_.size() * sizeof(T)));
        if (!out_) throw std::runtime_error("Scratch write failed (check free space): " + path_);
        buffer_.clear();
    }
    // ==========================================
    // Function: close
    // Method: Commit all records and release the write buffer.
    // ==========================================
    void close() {
        flush(); out_.close();
        if (!out_) throw std::runtime_error("Cannot close scratch file: " + path_);
        std::vector<T>().swap(buffer_);
    }
};

// ==========================================
// Class: BinaryReader
// Method: Read POD records in bounded blocks, rejecting truncated records.
// ==========================================
template<class T> class BinaryReader {
    std::string path_;
    std::ifstream in_;
    std::vector<T> buffer_;
    std::size_t position_ = 0, available_ = 0;
public:
    // ==========================================
    // Function: BinaryReader constructor
    // Method: Allocate only the configured read block.
    // ==========================================
    BinaryReader(const std::string &path, std::size_t bytes)
        : path_(path), in_(path, std::ios::binary),
          buffer_(std::max<std::size_t>(bytes / sizeof(T), 1)) {
        if (!in_) throw std::runtime_error("Cannot read scratch file: " + path);
    }
    // ==========================================
    // Function: next
    // Method: Return one record while reusing the fixed block allocation.
    // ==========================================
    bool next(T &value) {
        if (position_ == available_) {
            if (in_.eof()) return false;
            in_.read(reinterpret_cast<char *>(buffer_.data()),
                     static_cast<std::streamsize>(buffer_.size() * sizeof(T)));
            const auto bytes = in_.gcount();
            if (in_.bad() || (in_.fail() && !in_.eof()) || bytes % sizeof(T))
                throw std::runtime_error("Truncated or unreadable scratch file: " + path_);
            available_ = static_cast<std::size_t>(bytes) / sizeof(T); position_ = 0;
            if (!available_) return false;
        }
        value = buffer_[position_++]; return true;
    }
};

// ==========================================
// Class: BarRows
// Method: Supply repeated bounded scans to the unchanged stellar-bar fit logic.
// ==========================================
class BarRows {
    const std::vector<BarDatum> *memory_ = nullptr;
    std::string file_;
    int first_ = 0, last_ = 10000000;
public:
    // ==========================================
    // Function: BarRows constructors
    // Method: Select an existing in-memory table or a rank-local binary spool.
    // ==========================================
    explicit BarRows(const std::vector<BarDatum> &memory) : memory_(&memory) {}
    explicit BarRows(const std::string &file) : file_(file) {}
    // ==========================================
    // Function: subset
    // Method: Restrict a scan to a contiguous batch of physical exposure IDs.
    // ==========================================
    BarRows subset(int first, int last) const {
        auto result = *this; result.first_ = first; result.last_ = last; return result;
    }
    // ==========================================
    // Function: each
    // Method: Stream metadata without keeping all selected background rows in RAM.
    // ==========================================
    void each(const std::function<void(const BarDatum &)> &visit) const {
        if (memory_) {
            for (const auto &b : *memory_)
                if (b.exposure_id >= first_ && b.exposure_id <= last_) visit(b);
        } else {
            BinaryReader<BarDatum> input(file_, 1024 * 1024);
            BarDatum b;
            while (input.next(b))
                if (b.exposure_id >= first_ && b.exposure_id <= last_) visit(b);
        }
    }
};
} // namespace bcg
