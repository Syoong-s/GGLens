#pragma once

#include "bcg_catalog.hpp"

namespace bcg::archive {
// The durable formats use fixed-width integers and IEEE float64, with no struct padding.
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);

inline void hash_bytes(std::uint64_t &hash, const char *data, std::size_t size) {
    for (std::size_t i = 0; i < size; ++i) {
        hash ^= static_cast<unsigned char>(data[i]);
        hash *= 1099511628211ULL;
    }
}
inline constexpr std::uint64_t hash_start = 14695981039346656037ULL;

class Reader {
    std::ifstream input_;
    std::uintmax_t remaining_;
    std::uint64_t hash_ = hash_start;
public:
    explicit Reader(const fs::path &path) : input_(path, std::ios::binary | std::ios::ate) {
        if (!input_) throw std::runtime_error("Cannot open archive: " + path.string());
        const auto size = input_.tellg();
        if (size < 0 || !input_.seekg(0)) throw std::runtime_error("Cannot size archive: " + path.string());
        remaining_ = static_cast<std::uintmax_t>(size);
    }

    std::uintmax_t remaining() const { return remaining_; }
    std::uint64_t checksum() const { return hash_; }

    template<class T> void read(T *data, std::size_t count) {
        static_assert(std::is_arithmetic_v<T>);
        if (count > remaining_ / sizeof(T)) throw std::runtime_error("Truncated archive");
        if (!count) return;
        const auto bytes = count * sizeof(T);
        auto *raw = reinterpret_cast<char *>(data);
        if (!input_.read(raw, bytes)) throw std::runtime_error("Cannot read archive");
        remaining_ -= bytes;
        hash_bytes(hash_, raw, bytes);
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        for (std::size_t i = 0; i < count; ++i) std::reverse(raw + i*sizeof(T), raw + (i+1)*sizeof(T));
#endif
    }

    template<class T> T get() {
        T value;
        read(&value, 1);
        return value;
    }

    std::uint64_t verify() {
        const auto expected = checksum();
        if (get<std::uint64_t>() != expected || remaining())
            throw std::runtime_error("Invalid archive checksum or length");
        return expected;
    }
};

class Writer {
    std::ofstream output_;
    std::uint64_t hash_ = hash_start;
public:
    explicit Writer(const fs::path &path) : output_(path, std::ios::binary | std::ios::trunc) {
        if (!output_) throw std::runtime_error("Cannot write archive: " + path.string());
        output_.exceptions(std::ios::badbit | std::ios::failbit);
    }

    template<class T> void write(const T *data, std::size_t count) {
        static_assert(std::is_arithmetic_v<T>);
        if (!count) return;
        const auto *raw = reinterpret_cast<const char *>(data);
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        for (std::size_t i = 0; i < count; ++i) {
            std::array<char, sizeof(T)> bytes;
            std::reverse_copy(raw + i*sizeof(T), raw + (i+1)*sizeof(T), bytes.begin());
            output_.write(bytes.data(), bytes.size());
            hash_bytes(hash_, bytes.data(), bytes.size());
        }
#else
        output_.write(raw, count * sizeof(T));
        hash_bytes(hash_, raw, count * sizeof(T));
#endif
    }

    template<class T> void put(T value) { write(&value, 1); }

    void seal() {
        put(hash_);
        output_.close();
    }
};

inline bool has_magic(const fs::path &path, const std::string &magic) {
    std::ifstream input(path, std::ios::binary);
    std::string prefix(magic.size(), '\0');
    return input.read(prefix.data(), prefix.size()) && prefix == magic;
}

template<class Write> void publish(const fs::path &path, Write write) {
    const auto parent = fs::absolute(path).parent_path();
    fs::create_directories(parent);
    std::string temporary = (parent / ".fq-esd.XXXXXX").string();
    if (!mkdtemp(temporary.data())) throw std::runtime_error("Cannot stage archive in " + parent.string());
    const auto product = fs::path(temporary) / "product";
    std::error_code ignored;
    try {
        write(product);
        fs::rename(product, path);
    } catch (...) {
        fs::remove(product, ignored);
        fs::remove(temporary, ignored);
        throw;
    }
    fs::remove(temporary, ignored);
}
} // namespace bcg::archive
