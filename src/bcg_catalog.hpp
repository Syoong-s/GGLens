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
struct Source {
    double ra = 0, dec = 0, zp = 0;
    double g1 = 0, g2 = 0, de = 0, h1 = 0, h2 = 0;
    double chi = 0;
    double xyz[3] = {0, 0, 0};
    double gf1 = 0, gf2 = 0;
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
// Structure: Cuts
// Method: Apply shared FD quality thresholds without PT-column selection.
// ==========================================
struct Cuts {
    double snrfcut, snrlow, snrhigh, starcut;
    double chi2_thresh, flagcut, zplow, zphigh;
    double gf_lim, imaxcut, jmaxcut;
    double mag_low, mag_high, psf_chi2_mltp;
    double r_half_thresh, ft_cut, fg_cut, gold_cut, ext_cut;
    double star_bar_mltp, expo_bar_mltp;
    int max_dup;
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
// Function: read_group_manifest
// Method: Parse numbered (group, galaxy, member) blocks from one list file.
// ==========================================
inline std::vector<CatalogTriple> read_group_manifest(const std::string &path) {
    const fs::path list_path = canonical_catalog_file(path, "Group list");
    std::ifstream in(list_path);
    if (!in) throw std::runtime_error("Cannot open group list: " + path);
    const fs::path base = list_path.parent_path();
    std::string line;
    while (std::getline(in, line) && trim(line).empty()) {}
    if (trim(line) != "# group_cat galaxy_cat member_cat")
        throw std::runtime_error("Group list must start with '# group_cat galaxy_cat member_cat'");
    std::vector<std::string> fields;
    while (std::getline(in, line)) {
        line = trim(line);
        if (!line.empty() && line[0] != '#') fields.push_back(line);
    }
    if (!in.eof()) throw std::runtime_error("I/O failure reading group list: " + path);
    if (fields.empty() || fields.size() % 4 != 0)
        throw std::runtime_error("Group list needs numbered blocks of three catalog paths");
    std::unordered_set<int> ids;
    std::unordered_set<std::string> group_paths;
    std::vector<CatalogTriple> triples;
    triples.reserve(fields.size() / 4);
    for (std::size_t i = 0; i < fields.size(); i += 4) {
        std::size_t end = 0;
        long long id = 0;
        try { id = std::stoll(fields[i], &end, 10); }
        catch (const std::exception &) { throw std::runtime_error("Invalid group-list block number: " + fields[i]); }
        if (end != fields[i].size() || id <= 0 || id > std::numeric_limits<int>::max() ||
            !ids.insert(static_cast<int>(id)).second)
            throw std::runtime_error("Group-list block numbers must be unique positive integers");
        CatalogTriple triple;
        triple.id = static_cast<int>(id);
        triple.group = canonical_catalog_file(fields[i + 1], "Group catalog", base);
        triple.galaxy = canonical_catalog_file(fields[i + 2], "Galaxy catalog", base);
        triple.member = canonical_catalog_file(fields[i + 3], "Member catalog", base);
        if (!group_paths.insert(triple.group).second)
            throw std::runtime_error("Group catalog repeated in list: " + triple.group);
        triples.push_back(std::move(triple));
    }
    return triples;
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
// Function: catalog_fingerprint
// Method: Invalidate a selected-lens cache whenever input identity or cuts change.
// ==========================================
inline std::string catalog_fingerprint(const std::vector<CatalogTriple> &triples,
                                      double zmin, double zmax, double mmin, double mmax) {
    std::ostringstream out;
    out << std::setprecision(17) << "BCG_LENSES_V2|" << zmin << '|' << zmax
        << '|' << mmin << '|' << mmax;
    for (const auto &triple : triples) {
        out << "|block:" << triple.id;
        for (const auto &p : {triple.group, triple.galaxy, triple.member})
            out << '|' << p << ':' << fs::file_size(p) << ':'
                << fs::last_write_time(p).time_since_epoch().count();
    }
    return out.str();
}

// ==========================================
// Function: read_lens_cache
// Method: Load complete selected BCG lenses only when the fingerprint agrees.
// ==========================================
inline bool read_lens_cache(const std::string &path, const std::string &key,
                            const std::vector<CatalogTriple> &triples,
                            std::vector<Lens> &lenses) {
    if (path.empty() || !fs::exists(path)) return false;
    std::ifstream in(path);
    std::string stored;
    if (!std::getline(in, stored) || stored != key) return false;
    std::size_t n = 0;
    if (!(in >> n) || n == 0 || n > 100000000) return false;
    std::unordered_set<int> valid_ids;
    for (const auto &triple : triples) valid_ids.insert(triple.id);
    std::vector<Lens> temp(n);
    for (auto &v : temp) {
        if (!(in >> v.catalog_id >> v.id >> v.ra >> v.dec >> v.z >> v.log_mh)) return false;
        if (!valid_ids.count(v.catalog_id) || v.id <= 0 ||
            !std::isfinite(v.ra) || !std::isfinite(v.dec) ||
            !std::isfinite(v.z) || !std::isfinite(v.log_mh) ||
            v.ra < 0 || v.ra >= 360 || v.dec < -90 || v.dec > 90) return false;
    }
    lenses.swap(temp);
    return true;
}

// ==========================================
// Function: write_lens_cache
// Method: Write a complete cache to a process-unique temp file, then rename.
// ==========================================
inline void write_lens_cache(const std::string &path, const std::string &key,
                             const std::vector<Lens> &lenses) {
    if (path.empty()) return;
    const std::string tmp = path + ".tmp." + std::to_string(getpid());
    {
        std::ofstream out(tmp, std::ios::trunc);
        if (!out) throw std::runtime_error("Cannot write lens cache: " + tmp);
        out << key << '\n' << lenses.size() << '\n' << std::setprecision(17);
        for (const auto &v : lenses)
            out << v.catalog_id << ' ' << v.id << ' ' << v.ra << ' ' << v.dec
                << ' ' << v.z << ' ' << v.log_mh << '\n';
        out.flush();
        if (!out) throw std::runtime_error("Error writing lens cache: " + tmp);
    }
    fs::rename(tmp, path);
}

// ==========================================
// Function: select_group_lenses
// Method: Stream only needed columns and keep lenses passing z and mass cuts.
// ==========================================
inline std::vector<Lens> select_group_lenses(const std::vector<CatalogTriple> &triples,
                                              double zmin, double zmax,
                                              double mmin, double mmax) {
    std::vector<Lens> result;
    for (const auto &triple : triples) {
        std::ifstream in(triple.group);
        if (!in) throw std::runtime_error("Cannot read group catalog: " + triple.group);
        std::string line;
        while (std::getline(in, line)) {
            const char *p = line.c_str();
            char *end = nullptr;
            const auto id = std::strtoll(p, &end, 10); if (end == p) continue; p = end;
            std::strtod(p, &end); if (end == p) throw std::runtime_error("Bad group row: " + triple.group); p = end;
            std::strtod(p, &end); if (end == p) throw std::runtime_error("Bad group row: " + triple.group); p = end;
            std::strtod(p, &end); if (end == p) throw std::runtime_error("Bad group row: " + triple.group); p = end;
            const double z = std::strtod(p, &end); if (end == p) throw std::runtime_error("Bad group row: " + triple.group); p = end;
            const double mh = std::strtod(p, &end); if (end == p) throw std::runtime_error("Bad group row: " + triple.group);
            if (z >= zmin && z < zmax && mh >= mmin && mh < mmax)
                result.push_back({id, 0, 0, z, mh, triple.id});
        }
        if (!in.eof()) throw std::runtime_error("I/O failure reading: " + triple.group);
    }
    if (result.empty()) throw std::runtime_error("No groups pass z and mass cuts");
    return result;
}

// ==========================================
// Function: assign_bcg_coordinates
// Method: Scan each numbered block's member and galaxy file once.
// ==========================================
inline void assign_bcg_coordinates(std::vector<Lens> &lenses,
                                   const std::vector<CatalogTriple> &triples) {
    for (const auto &triple : triples) {
        std::unordered_map<std::int64_t, std::size_t> wanted;
        for (std::size_t i = 0; i < lenses.size(); ++i) {
            if (lenses[i].catalog_id != triple.id) continue;
            if (!wanted.emplace(lenses[i].id, i).second)
                throw std::runtime_error("Duplicate selected group ID in block " +
                                         std::to_string(triple.id) + ": " + std::to_string(lenses[i].id));
        }
        if (wanted.empty()) continue;
        std::ifstream min(triple.member);
        if (!min) throw std::runtime_error("Cannot read BCG member catalog: " + triple.member);
        std::unordered_map<std::int64_t, std::size_t> galaxy_to_lens;
        std::vector<unsigned char> found(lenses.size(), 0);
        std::string line;
        while (std::getline(min, line)) {
            const char *p = line.c_str(); char *end = nullptr;
            const auto galaxy_id = std::strtoll(p, &end, 10); if (end == p) continue; p = end;
            const auto group_id = std::strtoll(p, &end, 10); if (end == p) throw std::runtime_error("Bad member row"); p = end;
            const auto brightest = std::strtol(p, &end, 10); if (end == p) throw std::runtime_error("Bad member row");
            if (brightest != 1) continue;
            const auto it = wanted.find(group_id);
            if (it == wanted.end()) continue;
            if (found[it->second]++) throw std::runtime_error("Multiple BCGs for group " + std::to_string(group_id));
            if (galaxy_id <= 0 || !galaxy_to_lens.emplace(galaxy_id, it->second).second)
                throw std::runtime_error("Invalid or duplicate BCG galaxy ID");
        }
        if (!min.eof()) throw std::runtime_error("I/O failure reading: " + triple.member);
        for (const auto &entry : wanted)
            if (found[entry.second] != 1)
                throw std::runtime_error("Missing BCG for group " + std::to_string(entry.first));

        std::ifstream gin(triple.galaxy);
        if (!gin) throw std::runtime_error("Cannot read galaxy catalog: " + triple.galaxy);
        std::int64_t row_id = 0;
        std::size_t remaining = galaxy_to_lens.size();
        while (remaining && std::getline(gin, line)) {
            ++row_id;
            const auto it = galaxy_to_lens.find(row_id);
            if (it == galaxy_to_lens.end()) continue;
            const char *p = line.c_str(); char *end = nullptr;
            const double ra = std::strtod(p, &end); if (end == p) throw std::runtime_error("Bad BCG RA"); p = end;
            const double dec = std::strtod(p, &end); if (end == p) throw std::runtime_error("Bad BCG Dec");
            if (!std::isfinite(ra) || !std::isfinite(dec) || ra < 0 || ra >= 360 || dec < -90 || dec > 90)
                throw std::runtime_error("Invalid BCG coordinates at row " + std::to_string(row_id));
            lenses[it->second].ra = ra;
            lenses[it->second].dec = dec;
            --remaining;
        }
        if (remaining) throw std::runtime_error("BCG galaxy IDs beyond end of: " + triple.galaxy);
    }
}

// ==========================================
// Function: broadcast_lenses
// Method: Broadcast only selected and BCG-centered POD records from rank zero.
// ==========================================
inline void broadcast_lenses(std::vector<Lens> &lenses, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    int n = rank == 0 ? static_cast<int>(lenses.size()) : 0;
    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    if (n <= 0) throw std::runtime_error("No selected lenses to broadcast");
    if (rank != 0) lenses.resize(n);
    const auto bytes = static_cast<std::size_t>(n) * sizeof(Lens);
    if (bytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Lens broadcast exceeds MPI count limit");
    MPI_Bcast(lenses.data(), static_cast<int>(bytes), MPI_BYTE, 0, comm);
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

// ==========================================
// Structure: ExternalColumns
// Method: Set one-based positions in the external prefix for each used field.
// ==========================================
struct ExternalColumns {
    int ra, dec, zp, mag;
    int flags_footprint, flags_foreground, flags_gold, ext_mash;
};

// ==========================================
// Structure: ShearRow
// Method: Keep only the measurements used by source selection and pairing.
// ==========================================
struct ShearRow {
    enum Field {
        ra, dec, zp, mag, flags_footprint, flags_foreground, flags_gold, ext_mash,
        poly_chi2, sigma, nstar, imax, jmax, half_light_flux, half_light_area, flag,
        psf_fwhm, snr_f, gf1, gf2, g1, g2, de, h1, h2, chi2, expo_num, count
    };
    std::array<double, count> values{};
};

// ==========================================
// Structure: ShearLayout
// Method: Store the first catalog width and a row-field destination for each input column.
// ==========================================
struct ShearLayout {
    std::vector<int> destinations;
};

ShearLayout read_shear_layout(const std::string &first_path, const ExternalColumns &external,
                              const Cuts &cuts, MPI_Comm comm);
void for_each_local_source(const std::vector<std::string> &paths, const Cuts &cuts,
                           const ShearLayout &layout, MPI_Comm comm,
                           const std::function<void(const Source &, const BarDatum &)> &visit);
LocalSources read_local_sources(const std::vector<std::string> &paths, const Cuts &cuts,
                                const ShearLayout &layout, MPI_Comm comm);
} // namespace bcg
