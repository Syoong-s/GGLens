#include "bcg_dr9.hpp"

namespace bcg {
// ==========================================
// Function: read_group_manifest
// Method: Parse numbered (group, galaxy, member) blocks from one list file.
// ==========================================
std::vector<CatalogTriple> read_group_manifest(const std::string &path) {
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
// Function: catalog_fingerprint
// Method: Invalidate the cache when the center, cuts, input files or layout change.
// ==========================================
std::string catalog_fingerprint(const std::vector<CatalogTriple> &triples,
                                double zmin, double zmax, double mmin, double mmax, Center center,
                                const LensColumns &columns) {
    std::ostringstream out;
    out << std::setprecision(17) << "DR9_LENSES_V3|" << center_name(center) << '|' << zmin << '|' << zmax
        << '|' << mmin << '|' << mmax;
    const auto write_format = [&out](const char *name, const CatalogTextFormat &text) {
        out << '|' << name << ':' << static_cast<int>(text.delimiter) << ':' << text.header;
    };
    const auto &[group, galaxy, member] = columns;
    write_format("group", group.text);
    out << ':' << group.id << ':' << group.redshift << ':' << group.halo_mass
        << ':' << static_cast<int>(group.mass_format);
    write_format("galaxy", galaxy.text);
    out << ':' << galaxy.id << ':' << galaxy.ra << ':' << galaxy.dec << ':' << galaxy.luminosity;
    write_format("member", member.text);
    out << ':' << member.galaxy_id << ':' << member.group_id;
    if (center == Center::bcg) {
        out << "|bcg:" << galaxy.stellar_mass << ':' << static_cast<int>(galaxy.mass_format)
            << ':' << galaxy.photometry_select << ':' << member.brightest << ':' << member.brightest_value;
    }
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
// Method: Load complete selected lens centers only when the fingerprint agrees.
// ==========================================
bool read_lens_cache(const std::string &path, const std::string &key,
                     const std::vector<CatalogTriple> &triples, Center center,
                     std::vector<Lens> &lenses) {
    if (path.empty() || !fs::exists(path)) return false;
    std::ifstream in(path);
    std::string stored;
    if (!std::getline(in, stored) || stored != key) return false;
    std::size_t n = 0;
    if (!(in >> n) || n == 0 || n > std::numeric_limits<int>::max() / sizeof(Lens)) return false;
    std::unordered_set<int> valid_ids;
    for (const auto &triple : triples) valid_ids.insert(triple.id);
    std::vector<Lens> temp(n);
    for (auto &v : temp) {
        if (!(in >> v.catalog_id >> v.id >> v.ra >> v.dec >> v.z >> v.log_mh
                  >> v.galaxy_id >> v.stellar_mass >> v.photometry_select)) return false;
        if (!valid_ids.count(v.catalog_id) || v.id <= 0 ||
            !std::isfinite(v.ra) || !std::isfinite(v.dec) ||
            !std::isfinite(v.z) || !std::isfinite(v.log_mh) ||
            v.ra < 0 || v.ra >= 360 || v.dec < -90 || v.dec > 90) return false;
        if (center == Center::bcg) {
            if (v.galaxy_id < 0 || !(v.stellar_mass > 0 && std::isfinite(v.stellar_mass))) return false;
        } else if (v.galaxy_id != -1 || v.stellar_mass != 0 || v.photometry_select != 0) {
            return false;
        }
    }
    in >> std::ws;
    if (!in.eof() || in.bad()) return false;
    lenses.swap(temp);
    return true;
}

// ==========================================
// Function: write_lens_cache
// Method: Publish a complete cache using a unique temporary name on the shared filesystem.
// ==========================================
void write_lens_cache(const std::string &path, const std::string &key,
                      const std::vector<Lens> &lenses) {
    if (path.empty()) return;
    std::string staging = path + ".tmp.XXXXXX";
    if (!mkdtemp(staging.data()))
        throw std::runtime_error("Cannot create temporary lens cache directory for " + path);
    const auto tmp = fs::path(staging) / "cache";
    std::error_code ignored;
    try {
        std::ofstream out(tmp);
        if (!out) throw std::runtime_error("Cannot write lens cache: " + tmp.string());
        out << key << '\n' << lenses.size() << '\n' << std::setprecision(17);
        for (const auto &v : lenses)
            out << v.catalog_id << ' ' << v.id << ' ' << v.ra << ' ' << v.dec
                << ' ' << v.z << ' ' << v.log_mh << ' ' << v.galaxy_id
                << ' ' << v.stellar_mass << ' ' << v.photometry_select << '\n';
        out.close();
        if (!out) throw std::runtime_error("Error writing lens cache: " + tmp.string());
        fs::rename(tmp, path);
    } catch (...) {
        fs::remove(tmp, ignored);
        fs::remove(staging, ignored);
        throw;
    }
    fs::remove(staging, ignored);
}

} // namespace bcg
