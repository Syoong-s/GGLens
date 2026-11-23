#pragma once
#include "bcg_dr9.hpp"
#include "bcg_mpi_rows.hpp"
#include <charconv>
#include <string_view>

namespace bcg {
inline constexpr std::string_view catalog_whitespace = " \t\r";

inline bool is_catalog_data_row(std::string_view line) {
    const auto first = line.find_first_not_of(catalog_whitespace);
    return first != std::string_view::npos && line[first] != '#';
}

// Return the next field without copying. A null view marks the end of the row.
inline std::string_view take_catalog_field(std::string_view &text, CatalogDelimiter delimiter) {
    if (delimiter == CatalogDelimiter::whitespace) {
        const auto start = text.find_first_not_of(catalog_whitespace);
        if (start == std::string_view::npos) return {};
        text.remove_prefix(start);
        const auto field = text.substr(0, text.find_first_of(catalog_whitespace));
        text.remove_prefix(field.size());
        return field;
    }

    if (!text.data()) return {};
    std::size_t start = 0;
    if (delimiter == CatalogDelimiter::comma) {
        const auto first = text.find_first_not_of(catalog_whitespace);
        if (first != std::string_view::npos && text[first] == '"') {
            auto quote = text.find('"', first + 1);
            while (quote != std::string_view::npos && quote + 1 < text.size() && text[quote + 1] == '"')
                quote = text.find('"', quote + 2);
            start = quote == std::string_view::npos ? text.size() : quote + 1;
        }
    }
    const auto end = text.find(static_cast<char>(delimiter), start);
    const auto field = text.substr(0, end);
    text = end == std::string_view::npos ? std::string_view{} : text.substr(end + 1);
    return field;
}

inline std::string_view trim_catalog_field(std::string_view field) {
    const auto first = field.find_first_not_of(catalog_whitespace);
    if (first == std::string_view::npos) return {};
    return field.substr(first, field.find_last_not_of(catalog_whitespace) - first + 1);
}

// Index fields only as far as requested and reuse the span storage between rows.
class Dr9Fields {
    std::string_view remaining_;
    std::vector<std::string_view> fields_;
    const std::string &path_;
    CatalogDelimiter delimiter_;
public:
    Dr9Fields(const std::string &path, CatalogDelimiter delimiter) : path_(path), delimiter_(delimiter) {
        fields_.reserve(16);
    }

    void reset(std::string_view line) {
        remaining_ = line;
        fields_.clear();
    }

    template<class T> T number(int column) {
        while (fields_.size() < static_cast<std::size_t>(column)) {
            const auto field = take_catalog_field(remaining_, delimiter_);
            if (!field.data())
                throw std::runtime_error("Missing column " + std::to_string(column) + " in " + path_);
            fields_.push_back(field);
        }
        auto field = fields_[column - 1];
        if (delimiter_ != CatalogDelimiter::whitespace) {
            field = trim_catalog_field(field);
            if (delimiter_ == CatalogDelimiter::comma && field.size() >= 2 &&
                field.front() == '"' && field.back() == '"') {
                field = trim_catalog_field(field.substr(1, field.size() - 2));
            }
        }
        if (!field.empty()) {
            T value{};
            const auto end = field.data() + field.size();
            const auto start = field.data() + (field.front() == '+');
            const auto parsed = std::from_chars(start, end, value);
            if (parsed.ec == std::errc() && parsed.ptr == end) return value;
        }
        throw std::runtime_error("Invalid numeric value in " + path_ + ", column " +
                                 std::to_string(column) + ": " + std::string(field));
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
    std::int64_t row_number_ = 0;
public:
    // ==========================================
    // Function: Dr9Rows
    // Method: Open the range and advance past a partial leading line when necessary.
    // ==========================================
    Dr9Rows(const std::string &path, bool header, MPI_Comm comm) : path_(path) {
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
        if (!input_) throw std::runtime_error("Cannot open catalog: " + path);
        if (begin == 0) {
            if (header) {
                input_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                offset_ = input_.gcount();
            }
            return;
        }

        input_.seekg(static_cast<std::streamoff>(begin - 1));
        if (input_.get() == '\n') return;
        std::string partial;
        if (std::getline(input_, partial)) offset_ += partial.size() + 1;
        else if (input_.eof() && !input_.bad()) offset_ = end_;
        else throw std::runtime_error("Cannot align catalog input at " + path);
    }

    // ==========================================
    // Function: next
    // Method: Return the next data line and retain the offset for the following batch.
    // ==========================================
    bool next(std::string &line) {
        while (offset_ < end_ && std::getline(input_, line)) {
            offset_ += line.size() + 1;
            if (is_catalog_data_row(line)) {
                ++row_number_;
                return true;
            }
        }
        if (input_.bad() || (input_.fail() && !input_.eof()))
            throw std::runtime_error("Read failed for catalog: " + path_);
        return false;
    }

    std::int64_t row_number() const { return row_number_; }
};

std::int64_t preceding_catalog_rows(std::int64_t local_rows, MPI_Comm comm);
std::int64_t galaxy_row_offset(const std::string &path, const GalaxyColumns &columns, MPI_Comm comm);

// ==========================================
// Function: read_dr9_range
// Method: Stream the rank's catalog rows to a format-specific reader.
// ==========================================
template<class Visit>
std::int64_t read_dr9_range(const std::string &path, const CatalogTextFormat &format,
                             MPI_Comm comm, Visit visit) {
    Dr9Rows input(path, format.header, comm);
    Dr9Fields row(path, format.delimiter);
    std::string line;
    while (input.next(line)) {
        row.reset(line);
        visit(row, input.row_number());
    }
    return input.row_number();
}
} // namespace bcg
