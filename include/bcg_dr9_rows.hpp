#pragma once
#include "bcg_mpi_rows.hpp"
#include <charconv>
#include <string_view>

namespace bcg {
// ==========================================
// Class: Fields
// Method: Walk whitespace-delimited fields and parse only values the matching stage uses.
// ==========================================
class Dr9Fields {
    std::string_view text_;
    const std::string &path_;
public:
    Dr9Fields(const std::string &line, const std::string &path) : text_(line), path_(path) {}

    // ==========================================
    // Function: take
    // Method: Advance over one required field without allocating a string.
    // ==========================================
    std::string_view take() {
        constexpr std::string_view space = " \t\r";
        const auto start = text_.find_first_not_of(space);
        if (start == std::string_view::npos) throw std::runtime_error("Short DR9 row in " + path_);
        text_.remove_prefix(start);
        const auto end = std::min(text_.find_first_of(space), text_.size());
        const auto field = text_.substr(0, end);
        text_.remove_prefix(end);
        return field;
    }

    // ==========================================
    // Function: number
    // Method: Parse a complete ASCII numeric token and identify its catalog on error.
    // ==========================================
    template<class T> T number() {
        T value{};
        const auto field = take();
        auto end = field.data() + field.size();
        const auto start = field.data() + (field.front() == '+');
        const auto parsed = std::from_chars(start, end, value);
        if (parsed.ec != std::errc() || parsed.ptr != end)
            throw std::runtime_error("Invalid DR9 number in " + path_ + ": " + std::string(field));
        return value;
    }
};

// ==========================================
// Class: Dr9Rows
// Method: Give each MPI rank the lines that start in its assigned byte range.
// ==========================================
class Dr9Rows {
    std::vector<char> buffer_;
    std::ifstream input_;
    const std::string &path_;
    std::uint64_t offset_ = 0, end_ = 0;
public:
    // ==========================================
    // Function: Dr9Rows
    // Method: Open the range and advance past a partial leading line when necessary.
    // ==========================================
    Dr9Rows(const std::string &path, MPI_Comm comm) : path_(path) {
        int rank = 0, ranks = 1;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &ranks);
        const auto [begin, end] = partition(fs::file_size(path), rank, ranks);
        offset_ = begin;
        end_ = end;
        if (offset_ == end_) return;

        buffer_.resize(1 << 20);
        input_.rdbuf()->pubsetbuf(buffer_.data(), buffer_.size());
        input_.open(path, std::ios::binary);
        if (!input_) throw std::runtime_error("Cannot open DR9 catalog: " + path);
        if (begin == 0) return;

        input_.seekg(static_cast<std::streamoff>(begin - 1));
        if (input_.get() == '\n') return;
        std::string partial;
        if (std::getline(input_, partial)) offset_ += partial.size() + 1;
        else if (input_.eof() && !input_.bad()) offset_ = end_;
        else throw std::runtime_error("Cannot align DR9 input at " + path);
    }

    // ==========================================
    // Function: next
    // Method: Return the next data line and retain the offset for the following batch.
    // ==========================================
    bool next(std::string &line) {
        while (offset_ < end_ && std::getline(input_, line)) {
            offset_ += line.size() + 1;
            const auto first = line.find_first_not_of(" \t\r");
            if (first != std::string::npos && line[first] != '#') return true;
        }
        if (input_.bad() || (input_.fail() && !input_.eof()))
            throw std::runtime_error("Read failed for DR9 catalog: " + path_);
        return false;
    }
};

// ==========================================
// Function: read_dr9_range
// Method: Stream the rank's catalog rows to a format-specific reader.
// ==========================================
template<class Visit>
void read_dr9_range(const std::string &path, MPI_Comm comm, Visit visit) {
    Dr9Rows input(path, comm);
    std::string line;
    while (input.next(line)) visit(Dr9Fields(line, path));
}
} // namespace bcg
