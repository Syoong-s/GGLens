#include "bcg_archive.hpp"
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

namespace {
constexpr char cache_magic[] = "FQELENS1";

std::vector<Lens> read_binary_cache(const std::string &path, const std::string &key, Center center) {
    archive::Reader in(path);
    char prefix[sizeof(cache_magic) - 1];
    in.read(prefix, sizeof(prefix));
    if (std::string(prefix, sizeof(prefix)) != cache_magic || in.get<std::uint64_t>() != key.size()) return {};
    std::string stored(key.size(), '\0');
    in.read(stored.data(), stored.size());
    if (stored != key) return {};
    const auto n = in.get<std::uint32_t>();
    const auto record_size = center == Center::bcg ? 64ULL : 44ULL;
    if (n == 0 || n > std::numeric_limits<int>::max() / sizeof(Lens) ||
        in.remaining() != record_size*n + sizeof(std::uint64_t)) return {};
    std::vector<Lens> lenses(n);
    for (auto &lens : lenses) {
        lens.catalog_id = in.get<std::int32_t>();
        lens.id = in.get<std::int64_t>();
        lens.ra = in.get<double>();
        lens.dec = in.get<double>();
        lens.z = in.get<double>();
        lens.log_mh = in.get<double>();
        if (center == Center::bcg) {
            lens.galaxy_id = in.get<std::int64_t>();
            lens.stellar_mass = in.get<double>();
            lens.photometry_select = in.get<std::int32_t>();
        }
    }
    in.verify();
    return lenses;
}

std::vector<Lens> read_text_cache(const std::string &path, const std::string &key) {
    std::ifstream in(path);
    std::string stored;
    if (!std::getline(in, stored) || stored != key) return {};
    std::size_t n = 0;
    if (!(in >> n) || n == 0 || n > std::numeric_limits<int>::max() / sizeof(Lens)) return {};
    std::vector<Lens> lenses(n);
    for (auto &lens : lenses)
        if (!(in >> lens.catalog_id >> lens.id >> lens.ra >> lens.dec >> lens.z >> lens.log_mh
                 >> lens.galaxy_id >> lens.stellar_mass >> lens.photometry_select)) return {};
    in >> std::ws;
    if (!in.eof() || in.bad()) return {};
    return lenses;
}
} // namespace

bool read_lens_cache(const std::string &path, const std::string &key,
                     const std::vector<CatalogTriple> &triples, Center center,
                     std::vector<Lens> &lenses) {
    if (path.empty() || !fs::is_regular_file(path)) return false;
    std::vector<Lens> cached;
    try {
        cached = archive::has_magic(path, cache_magic) ? read_binary_cache(path, key, center) :
                                                       read_text_cache(path, key);
    } catch (const std::runtime_error &) {
        return false;
    }
    if (cached.empty()) return false;
    std::unordered_set<int> valid_ids;
    for (const auto &triple : triples) valid_ids.insert(triple.id);
    for (const auto &lens : cached) {
        if (!valid_ids.count(lens.catalog_id) || lens.id <= 0 ||
            !std::isfinite(lens.ra) || !std::isfinite(lens.dec) ||
            !std::isfinite(lens.z) || !std::isfinite(lens.log_mh) ||
            lens.ra < 0 || lens.ra >= 360 || lens.dec < -90 || lens.dec > 90) return false;
        if (center == Center::bcg) {
            if (lens.galaxy_id < 0 || !(lens.stellar_mass > 0 && std::isfinite(lens.stellar_mass))) return false;
        } else if (lens.galaxy_id != -1 || lens.stellar_mass != 0 || lens.photometry_select != 0) {
            return false;
        }
    }
    lenses.swap(cached);
    return true;
}

void write_lens_cache(const std::string &path, const std::string &key, Center center,
                      const std::vector<Lens> &lenses) {
    if (path.empty()) return;
    archive::publish(path, [&](const auto &temporary) {
        archive::Writer out(temporary);
        out.write(cache_magic, sizeof(cache_magic) - 1);
        out.put(static_cast<std::uint64_t>(key.size()));
        out.write(key.data(), key.size());
        out.put(static_cast<std::uint32_t>(lenses.size()));
        for (const auto &lens : lenses) {
            out.put(static_cast<std::int32_t>(lens.catalog_id));
            out.put(lens.id);
            out.put(lens.ra);
            out.put(lens.dec);
            out.put(lens.z);
            out.put(lens.log_mh);
            if (center == Center::bcg) {
                out.put(lens.galaxy_id);
                out.put(lens.stellar_mass);
                out.put(static_cast<std::int32_t>(lens.photometry_select));
            }
        }
        out.seal();
    });
}
} // namespace bcg
