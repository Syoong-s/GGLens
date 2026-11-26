#pragma once

#include <mpi.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unistd.h>

namespace bcg {
namespace fs = std::filesystem;

struct Lens {
    std::int64_t id = 0;
    double ra = 0, dec = 0, z = 0, log_mh = 0;
    int catalog_id = 0; // Numbered block in the group input list.
    std::int64_t galaxy_id = -1;
    double stellar_mass = 0; // Linear physical M_sun.
    int photometry_select = 0;
};
static_assert(std::is_trivially_copyable<Lens>::value, "MPI Lens broadcast requires POD");

struct CatalogTriple {
    int id = 0;
    std::string group, galaxy, member;
};

// ==========================================
// Structure: Source
// Method: Retain Cartesian moments and field distortion for source-level calibration.
// ==========================================
struct SourceRecord {
    double ra = 0, dec = 0, zp = 0;
    double g1 = 0, g2 = 0, de = 0, h1 = 0, h2 = 0;
    double gf1 = 0, gf2 = 0;
    double zp_err = 0;
};

struct Source : SourceRecord {
    double chi = 0;
    double xyz[3] = {0, 0, 0};

    Source() = default;
    explicit Source(const SourceRecord &record) : SourceRecord(record) {}
};

// ==========================================
// Structure: BarDatum
// Method: Retain only the fields needed for an optional stellar-bar fit.
// ==========================================
struct BarDatum {
    double mag_i = 0, size_rel = 0, snr = 0;
    int exposure_id = 0;
};

// ==========================================
// Structure: LocalSources
// Method: Keep bar metadata aligned with quality-selected source rows.
// ==========================================
struct LocalSources {
    std::vector<Source> sources;
    std::vector<BarDatum> bar_data;
};

// ==========================================
// Function: trim
// Method: Remove surrounding whitespace from a catalog-list entry.
// ==========================================
inline std::string trim(const std::string &s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// ==========================================
// Function: canonical_catalog_file
// Method: Resolve a list entry against its list directory and require a regular file.
// ==========================================
inline std::string canonical_catalog_file(const std::string &path, const std::string &role,
                                          const fs::path &base = {}) {
    std::string value = trim(path);
    if (value.size() >= 2 &&
        ((value.front() == '"' && value.back() == '"') ||
         (value.front() == '\'' && value.back() == '\'')))
        value = value.substr(1, value.size() - 2);
    fs::path candidate(value);
    if (candidate.is_relative()) candidate = base / candidate;
    if (!fs::is_regular_file(candidate))
        throw std::runtime_error(role + " is not a regular file: " + candidate.string());
    return fs::canonical(candidate).string();
}

// ==========================================
// Function: read_shear_list
// Method: Require a list file with one shear-catalog path per line.
// ==========================================
inline std::vector<std::string> read_shear_list(const std::string &path) {
    const fs::path list_path = canonical_catalog_file(path, "Shear list");
    std::ifstream in(list_path);
    if (!in) throw std::runtime_error("Cannot open shear list: " + path);
    const fs::path base = list_path.parent_path();
    std::vector<std::string> files;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        files.push_back(canonical_catalog_file(line, "Shear catalog", base));
    }
    if (!in.eof()) throw std::runtime_error("I/O failure reading shear list: " + path);
    if (files.empty()) throw std::runtime_error("Shear list contains no catalog paths: " + path);
    return files;
}

// ==========================================
// Function: broadcast_paths
// Method: Send one newline-delimited shear-file list in one collective.
// ==========================================
inline void broadcast_paths(std::vector<std::string> &paths, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    std::string packed;
    if (rank == 0) for (const auto &p : paths) packed += p + '\n';
    std::uint64_t length = packed.size();
    MPI_Bcast(&length, 1, MPI_UINT64_T, 0, comm);
    if (length > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Shear list exceeds MPI count limit");
    if (rank != 0) packed.resize(length);
    MPI_Bcast(packed.data(), static_cast<int>(length), MPI_CHAR, 0, comm);
    if (rank != 0) {
        std::istringstream in(packed);
        std::string p;
        while (std::getline(in, p)) if (!p.empty()) paths.push_back(p);
    }
}

} // namespace bcg
