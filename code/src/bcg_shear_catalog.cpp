#include "bcg_shear_cuts.hpp"
#include <string_view>

namespace bcg {
namespace {
using Field = ShearRow::Field;

// ==========================================
// Function: detect_shear_layout
// Method: Map the first header's shear.dat fields and the configured metadata positions.
// ==========================================
ShearLayout detect_shear_layout(const std::string &path, const ExternalColumns &external,
                                const Cuts &cuts, bool use_pair_redshift_error) {
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
    if (names.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) / sizeof(ShearLayout::Column))
        throw std::runtime_error("Shear header exceeds MPI count limit: " + path);

    ShearLayout layout;
    layout.columns.resize(names.size());
    const auto assign = [&](std::size_t index, Field field) {
        auto &column = layout.columns[index];
        if (column.destination >= 0)
            throw std::runtime_error("Shear column assigned twice: " + std::to_string(index + 1));
        column.destination = field;
    };
    const int width = static_cast<int>(names.size());
    const auto external_index = [&](int column, const char *key) {
        if (column < 1 || column > width)
            throw std::runtime_error("Set external_columns." + std::string(key) +
                " to a valid one-based column number for " + path);
        return static_cast<std::size_t>(column - 1);
    };
    const auto set_external = [&](int column, Field field, const char *key, bool required) {
        if (required) assign(external_index(column, key), field);
    };
    const auto set_limits = [&](int column, std::optional<double> low, std::optional<double> high, const char *key) {
        if (!low && !high) return;
        auto &mapping = layout.columns[external_index(column, key)];
        if (low) mapping.low = std::max(mapping.low, *low);
        if (high) mapping.high = std::min(mapping.high, *high);
        mapping.bounded = true;
    };
    set_external(external.ra, Field::ra, "ra", true);
    set_external(external.dec, Field::dec, "dec", true);
    set_external(external.zp, Field::zp, "zp", true);
    if (use_pair_redshift_error) {
        constexpr char key[] = "ext-zp-err";
        const auto rule = std::find_if(cuts.external.begin(), cuts.external.end(),
            [&](const auto &cut) { return cut.key == key; });
        const int column = rule == cuts.external.end() ? 0 : rule->column;
        set_external(column, Field::zp_err, key, true);
    }
    set_external(external.mag_stellar_locus, Field::mag_stellar_locus, "mag-stellar-locus", cuts.bar_enabled());
    set_external(external.expo_num, Field::expo_num, "expo-num", cuts.exposure_bar_enabled());
    set_external(external.ccd_num, Field::ccd_num, "ccd-num", cuts.ccd_cut_enabled());
    const int chi2_column = external.expo_chi2 == -1 ? width : external.expo_chi2;
    set_external(chi2_column, Field::chi2, "expo-chi2", cut_enabled(cuts.chi2_thresh));

    struct NamedColumn {
        std::string_view name;
        Field field;
        bool needed = true;
    };
    const NamedColumn shear_fields[] = {
        // Moments and field distortion are used by measurement.
        {"g1", Field::g1}, {"g2", Field::g2}, {"de", Field::de},
        {"h1", Field::h1}, {"h2", Field::h2}, {"gf1", Field::gf1}, {"gf2", Field::gf2},
        // Read quality columns when their cuts or the stellar bar consume them.
        {"snr_f", Field::snr_f, cut_enabled(cuts.snrfcut)},
        {"half_light_flux", Field::half_light_flux, cuts.needs_snr()},
        {"half_light_area", Field::half_light_area, cuts.needs_snr() || cuts.needs_psf()},
        {"psf_fwhm", Field::psf_fwhm, cuts.needs_psf()},
        {"nstar", Field::nstar, cut_enabled(cuts.starcut)},
        {"poly_chi2", Field::poly_chi2, cut_enabled(cuts.psf_chi2_mltp)},
        {"flag", Field::flag, cut_enabled(cuts.flagcut)},
        {"sigma", Field::sigma, cut_enabled(cuts.flagcut)},
        {"imax", Field::imax, cut_enabled(cuts.imaxcut)},
        {"jmax", Field::jmax, cut_enabled(cuts.jmaxcut)},
    };
    // CatalogCombiner appends shear.dat after the external catalog. A repeated
    // label in the external data therefore precedes the shear.dat occurrence.
    for (const auto &[name, field, needed] : shear_fields) {
        if (!needed) continue;
        auto column = names.size();
        while (column > 0 && (names[column - 1] != name || layout.columns[column - 1].used()))
            --column;
        if (column == 0)
            throw std::runtime_error("Missing shear.dat field " + std::string(name) + ": " + path);
        assign(column - 1, field);
    }
    // Custom bounds can also use columns already mapped to named shear fields.
    for (const auto &cut : cuts.external) set_limits(cut.column, cut.low, cut.high, cut.key.c_str());
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
    for (std::size_t i = 0; i < layout.columns.size(); ++i) {
        const char *start = p;
        while (*p && *p != '#' && !std::isspace(static_cast<unsigned char>(*p))) ++p;
        const auto &column = layout.columns[i];
        if (column.used()) {
            char *end = nullptr;
            const double value = std::strtod(start, &end);
            if (end == start || end != p)
                throw std::runtime_error("Invalid numeric shear value at " + path + ":" +
                    std::to_string(lineno) + ", column " + std::to_string(i + 1));
            if (!std::isfinite(value) || std::abs(value) > 1e300) finite_row = false;
            if (column.destination >= 0) row.values[static_cast<std::size_t>(column.destination)] = value;
            // Leave finite rows in the stream so they still count toward duplicate groups.
            if (column.bounded && (value < column.low || value > column.high))
                row.passes_column_cuts = false;
        }
        while (std::isspace(static_cast<unsigned char>(*p))) ++p;
    }
    return finite_row;
}

} // namespace

// ==========================================
// Function: read_shear_layout
// Method: Read only the first catalog header on rank zero and broadcast its width and mapping.
// ==========================================
ShearLayout read_shear_layout(const std::string &first_path, const ExternalColumns &external,
                              const Cuts &cuts, bool use_pair_redshift_error, MPI_Comm comm) {
    static_assert(std::is_trivially_copyable_v<ShearLayout::Column>, "MPI shear mapping requires POD");
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    ShearLayout layout;
    if (rank == 0) layout = detect_shear_layout(first_path, external, cuts, use_pair_redshift_error);
    int width = static_cast<int>(layout.columns.size());
    MPI_Bcast(&width, 1, MPI_INT, 0, comm);
    if (rank != 0) layout.columns.resize(static_cast<std::size_t>(width));
    MPI_Bcast(layout.columns.data(), static_cast<int>(layout.columns.size() * sizeof(ShearLayout::Column)),
              MPI_BYTE, 0, comm);
    return layout;
}

// ==========================================
// Function: for_each_shear_row
// Method: Decode one catalog with the shared column mapping and emit its numeric rows in order.
// ==========================================
void for_each_shear_row(const std::string &path, const ShearLayout &layout,
                        const std::function<void(const ShearRow &)> &visit) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read shear catalog: " + path);
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::string line;
    std::size_t lineno = 1;
    while (std::getline(in, line)) {
        ++lineno;
        ShearRow row;
        if (parse_shear_values(line, row, layout, path, lineno)) visit(row);
    }
    if (!in.eof()) throw std::runtime_error("I/O failure reading: " + path);
}
} // namespace bcg
