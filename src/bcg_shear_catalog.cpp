#include "bcg_catalog.hpp"
#include <string_view>

namespace bcg {
namespace {
using Field = ShearRow::Field;

// ==========================================
// Function: detect_shear_layout
// Method: Count the first header's columns and resolve only the used FQ suffix fields.
// ==========================================
ShearLayout detect_shear_layout(const std::string &path, const ExternalColumns &external,
                                const Cuts &cuts) {
    std::ifstream in(path);
    std::string header;
    if (!std::getline(in, header)) throw std::runtime_error("Cannot read shear header: " + path);
    header = trim(header);
    if (!header.empty() && header.front() == '#') header.erase(0, 1);
    std::istringstream input(header);
    std::vector<std::string> names;
    for (std::string name; input >> name;) {
        std::transform(name.begin(), name.end(), name.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        names.push_back(std::move(name));
    }
    if (names.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Shear header exceeds MPI count limit: " + path);

    const auto find_column = [&](std::string_view name, auto begin) {
        const auto it = std::find(begin, names.cend(), name);
        if (it != names.cend() && std::find(it + 1, names.cend(), name) != names.cend())
            throw std::runtime_error("Ambiguous shear field " + std::string(name) + ": " + path);
        return it;
    };
    const auto exposure = find_column("expo_num", names.cbegin());
    const auto marker = exposure != names.cend() ? exposure : find_column("ccd_num", names.cbegin());
    if (marker == names.cend())
        throw std::runtime_error("Shear header requires EXPO_NUM or the legacy CCD_NUM: " + path);
    if (cuts.expo_bar_mltp > 0 && exposure == names.cend())
        throw std::runtime_error("Exposure bar requires EXPO_NUM: " + path);

    ShearLayout layout;
    layout.destinations.assign(names.size(), -1);
    const auto assign = [&](std::size_t index, Field field) {
        if (layout.destinations[index] != -1)
            throw std::runtime_error("Shear column assigned twice: " + std::to_string(index + 1));
        layout.destinations[index] = field;
    };
    const auto prefix_width = marker - names.cbegin();
    const auto set_external = [&](int column, Field field, const char *key, bool required) {
        if (!required) return;
        if (column < 1 || column > prefix_width)
            throw std::runtime_error("external_columns." + std::string(key) +
                " must point to the external prefix of " + path);
        assign(static_cast<std::size_t>(column - 1), field);
    };

    set_external(external.ra, Field::ra, "ra", true);
    set_external(external.dec, Field::dec, "dec", true);
    set_external(external.zp, Field::zp, "zp", true);
    set_external(external.mag, Field::mag, "mag", true);
    set_external(external.flags_footprint, Field::flags_footprint, "flags-footprint", cuts.ft_cut >= 0);
    set_external(external.flags_foreground, Field::flags_foreground, "flags-foreground", cuts.fg_cut >= 0);
    set_external(external.flags_gold, Field::flags_gold, "flags-gold", cuts.gold_cut >= 0);
    set_external(external.ext_mash, Field::ext_mash, "ext-mash", cuts.ext_cut >= 0);

    const std::pair<std::string_view, Field> required[] = {
        {"poly_chi2", Field::poly_chi2}, {"sigma", Field::sigma}, {"nstar", Field::nstar},
        {"imax", Field::imax}, {"jmax", Field::jmax},
        {"half_light_flux", Field::half_light_flux}, {"half_light_area", Field::half_light_area},
        {"flag", Field::flag}, {"psf_fwhm", Field::psf_fwhm}, {"snr_f", Field::snr_f},
        {"gf1", Field::gf1}, {"gf2", Field::gf2}, {"g1", Field::g1}, {"g2", Field::g2},
        {"de", Field::de}, {"h1", Field::h1}, {"h2", Field::h2}, {"chi2", Field::chi2}
    };
    for (const auto &[name, field] : required) {
        const auto column = find_column(name, marker + 1);
        if (column == names.cend())
            throw std::runtime_error("Missing FQ suffix field " + std::string(name) + ": " + path);
        assign(static_cast<std::size_t>(column - names.cbegin()), field);
    }
    if (exposure != names.cend())
        assign(static_cast<std::size_t>(exposure - names.cbegin()), Field::expo_num);
    return layout;
}

// ==========================================
// Function: parse_shear_values
// Method: Extract mapped numeric fields from the shared layout and skip unused tokens.
// ==========================================
bool parse_shear_values(const std::string &line, ShearRow &row, const ShearLayout &layout,
                        const std::string &path, std::size_t lineno) {
    const char *p = line.c_str();
    while (std::isspace(static_cast<unsigned char>(*p))) ++p;
    if (!*p || *p == '#') return false;
    bool finite_row = true;
    for (std::size_t i = 0; i < layout.destinations.size(); ++i) {
        const char *start = p;
        while (*p && *p != '#' && !std::isspace(static_cast<unsigned char>(*p))) ++p;
        const int destination = layout.destinations[i];
        if (destination >= 0) {
            char *end = nullptr;
            const double value = std::strtod(start, &end);
            if (end == start || end != p)
                throw std::runtime_error("Invalid numeric shear value at " + path + ":" +
                    std::to_string(lineno) + ", column " + std::to_string(i + 1));
            if (!std::isfinite(value) || std::abs(value) > 1e300) finite_row = false;
            row.values[static_cast<std::size_t>(destination)] = value;
        }
        while (std::isspace(static_cast<unsigned char>(*p))) ++p;
    }
    return finite_row;
}

// ==========================================
// Function: passes_shear_cuts
// Method: Apply FD quality cuts using mapped fields and retain celestial moments.
// ==========================================
bool passes_shear_cuts(const ShearRow &row, const Cuts &cuts, Source &s, BarDatum &bar) {
    const auto &v = row.values;
    if (cuts.ft_cut >= 0 && std::abs(v[Field::flags_footprint] - cuts.ft_cut) > 1e-3) return false;
    if (cuts.fg_cut >= 0 && std::abs(v[Field::flags_foreground] - cuts.fg_cut) > 1e-3) return false;
    if (cuts.gold_cut >= 0 && std::abs(v[Field::flags_gold] - cuts.gold_cut) > 1e-3) return false;
    if (cuts.ext_cut >= 0 && std::abs(v[Field::ext_mash] - cuts.ext_cut) > 1e-3) return false;
    const double area = v[Field::half_light_area], psf = v[Field::psf_fwhm];
    if (area <= 0 || psf <= 0 || v[Field::snr_f] < cuts.snrfcut) return false;
    const double snr = v[Field::half_light_flux] / std::sqrt(area);
    if (cuts.snrlow > 0 && snr < cuts.snrlow) return false;
    if (cuts.snrhigh > 0 && snr > cuts.snrhigh) return false;
    constexpr double pi = 3.14159265358979323846;
    if (cuts.r_half_thresh > 0) {
        const double rhalf = std::sqrt(area / pi) * 0.2628 * 2.0;
        if (rhalf <= cuts.r_half_thresh * psf) return false;
    }
    if (v[Field::mag] < cuts.mag_low || v[Field::mag] > cuts.mag_high) return false;
    if (v[Field::poly_chi2] > cuts.psf_chi2_mltp || v[Field::nstar] < cuts.starcut ||
        v[Field::chi2] > cuts.chi2_thresh) return false;
    if (cuts.flagcut > 0 && v[Field::flag] * v[Field::sigma] > cuts.flagcut) return false;
    if (v[Field::imax] >= cuts.imaxcut || v[Field::jmax] >= cuts.jmaxcut) return false;
    if (v[Field::zp] <= cuts.zplow || v[Field::zp] >= cuts.zphigh) return false;
    if (std::abs(v[Field::gf1]) > cuts.gf_lim || std::abs(v[Field::gf2]) > cuts.gf_lim) return false;
    if (v[Field::ra] < 0 || v[Field::ra] >= 360 || v[Field::dec] < -90 || v[Field::dec] > 90) return false;
    s = {v[Field::ra], v[Field::dec], v[Field::zp], v[Field::g1], v[Field::g2],
         v[Field::de], v[Field::h1], v[Field::h2]};
    s.gf1 = v[Field::gf1];
    s.gf2 = v[Field::gf2];
    if (cuts.star_bar_mltp >= 0 || cuts.expo_bar_mltp > 0) {
        bar.mag_i = v[Field::mag];
        bar.snr = snr;
        bar.size_rel = (2.0 * 0.2628 * std::sqrt(area / pi) - psf) / psf;
        if (cuts.expo_bar_mltp > 0) {
            const double exposure = v[Field::expo_num];
            if (exposure < 1 || exposure > 10000000)
                throw std::runtime_error("Invalid exposure ID for enabled exposure bar");
            bar.exposure_id = static_cast<int>(std::llround(exposure));
            if (std::abs(exposure - bar.exposure_id) > 0.01)
                throw std::runtime_error("Non-integer exposure ID for enabled exposure bar");
        }
    }
    return true;
}
} // namespace

// ==========================================
// Function: read_shear_layout
// Method: Read only the first catalog header on rank zero and broadcast its width and mapping.
// ==========================================
ShearLayout read_shear_layout(const std::string &first_path, const ExternalColumns &external,
                              const Cuts &cuts, MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    ShearLayout layout;
    if (rank == 0) layout = detect_shear_layout(first_path, external, cuts);
    int width = static_cast<int>(layout.destinations.size());
    MPI_Bcast(&width, 1, MPI_INT, 0, comm);
    if (rank != 0) layout.destinations.resize(static_cast<std::size_t>(width));
    MPI_Bcast(layout.destinations.data(), width, MPI_INT, 0, comm);
    return layout;
}

// ==========================================
// Function: for_each_local_source
// Method: Stream adjacent-position groups using the broadcast layout and skip all catalog headers.
// ==========================================
void for_each_local_source(const std::vector<std::string> &paths, const Cuts &cuts,
                           const ShearLayout &layout, MPI_Comm comm,
                           const std::function<void(const Source &, const BarDatum &)> &visit) {
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    for (std::size_t i = 0; i < paths.size(); ++i) {
        if (size > 1 && (rank == 0 || static_cast<int>(i % (size - 1)) != rank - 1)) continue;
        std::ifstream in(paths[i]);
        if (!in) throw std::runtime_error("Cannot read shear catalog: " + paths[i]);
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::string line;
        std::size_t lineno = 1, duplicate_count = 0;
        double last_ra = 0, last_dec = 0;
        std::vector<ShearRow> duplicates;
        duplicates.reserve(std::min<std::size_t>(static_cast<std::size_t>(cuts.max_dup), 64));
        const auto flush = [&] {
            if (duplicate_count <= static_cast<std::size_t>(cuts.max_dup))
                for (const auto &row : duplicates) {
                    Source source;
                    BarDatum bar;
                    if (passes_shear_cuts(row, cuts, source, bar)) visit(source, bar);
                }
            duplicates.clear();
            duplicate_count = 0;
        };
        while (std::getline(in, line)) {
            ++lineno;
            ShearRow row;
            if (!parse_shear_values(line, row, layout, paths[i], lineno)) continue;
            const double ra = row.values[Field::ra], dec = row.values[Field::dec];
            if (duplicate_count && (std::abs(ra - last_ra) > 1e-7 ||
                                    std::abs(dec - last_dec) > 1e-7)) flush();
            if (!duplicate_count) {
                last_ra = ra;
                last_dec = dec;
            }
            ++duplicate_count;
            if (duplicate_count <= static_cast<std::size_t>(cuts.max_dup))
                duplicates.push_back(row);
        }
        flush();
        if (!in.eof()) throw std::runtime_error("I/O failure reading: " + paths[i]);
    }
}

// ==========================================
// Function: read_local_sources
// Method: Collect selected rows and aligned optional bar metadata through the streaming reader.
// ==========================================
LocalSources read_local_sources(const std::vector<std::string> &paths, const Cuts &cuts,
                                const ShearLayout &layout, MPI_Comm comm) {
    LocalSources result;
    const bool need_bar = cuts.star_bar_mltp >= 0 || cuts.expo_bar_mltp > 0;
    for_each_local_source(paths, cuts, layout, comm, [&](const Source &s, const BarDatum &b) {
        result.sources.push_back(s);
        if (need_bar) result.bar_data.push_back(b);
    });
    return result;
}
} // namespace bcg
