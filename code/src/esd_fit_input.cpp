#include "bcg_esd.hpp"
#include "bcg_archive.hpp"
#include "selfcalibration.hpp"
#include <iostream>

namespace bcg_esd {
namespace {
constexpr char format_name[] = "FQ_GGLENS_FIT_INPUT";
constexpr int format_version = 1;

void expect(std::istream &in, const char *key) {
    std::string token;
    if (!(in >> token) || token != key)
        throw std::runtime_error(std::string("Expected ") + key);
}

template<class T>
void read_field(std::istream &in, const char *key, T &value) {
    expect(in, key);
    if (!(in >> value)) throw std::runtime_error(std::string("Invalid value for ") + key);
}

void expect_end(std::istream &in) {
    in >> std::ws;
    if (!in.eof() || in.bad()) throw std::runtime_error("Unexpected data after the last field");
}

std::ifstream open_input(const bcg::fs::path &path) {
    if (!bcg::fs::is_regular_file(path))
        throw std::runtime_error("ESD fit input is not a regular file: " + path.string());
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open ESD fit input: " + path.string());
    return in;
}

void validate_bins(const std::vector<int> &bins, int radial_bins) {
    int previous = -1;
    for (const int bin : bins) {
        if (bin <= previous || bin >= radial_bins)
            throw std::runtime_error("Fit bins must be distinct, ordered indices within the measurement");
        previous = bin;
    }
}

FitInput read_snapshot(std::istream &in) {
    int version = 0;
    read_field(in, format_name, version);
    if (version != format_version) throw std::runtime_error("Unsupported ESD fit input version");

    FitInput input;
    auto &[center, cosmology, sample, radial_bins] = input.metadata;
    std::string name;
    read_field(in, "center", name);
    center = bcg::parse_center(name);
    expect(in, "sample");
    if (!(in >> sample.n_lenses >> sample.regions >> sample.redshift >> sample.stellar_mass))
        throw std::runtime_error("Invalid lens sample metadata");
    expect(in, "cosmology");
    if (!(in >> cosmology.omega_m >> cosmology.h >> cosmology.sigma8 >> cosmology.omega_b_h2
             >> cosmology.spectral_index >> cosmology.neutrino_mass))
        throw std::runtime_error("Invalid cosmology metadata");
    int count = 0;
    read_field(in, "radial_bins", radial_bins);
    read_field(in, "fit_bins", count);
    if (radial_bins < 1 || radial_bins > max_radial_bins || count < 0 || count > radial_bins)
        throw std::runtime_error("Invalid radial and fitted bin counts");

    input.bins.resize(count);
    input.radius.resize(count);
    input.measured.resize(count);
    input.covariance.resize(count * count);
    expect(in, "index_radius_esd");
    for (int i = 0; i < count; ++i) {
        if (!(in >> input.bins[i] >> input.radius[i] >> input.measured[i]))
            throw std::runtime_error("Invalid ESD measurement row");
    }
    expect(in, "covariance");
    for (double &value : input.covariance) {
        if (!(in >> value)) throw std::runtime_error("Invalid fitted covariance matrix");
    }
    expect_end(in);
    return input;
}

std::istringstream legacy_header(std::string line) {
    std::replace(line.begin(), line.end(), '=', ' ');
    return std::istringstream(line.substr(2));
}

// Saved tables use nan for unavailable values. Formatted double extraction does not read it.
double table_number(std::istream &in) {
    std::string token;
    if (!(in >> token)) throw std::runtime_error("Incomplete table row");
    char *end = nullptr;
    const double value = std::strtod(token.c_str(), &end);
    if (end != token.c_str() + token.size()) throw std::runtime_error("Invalid table value: " + token);
    return value;
}

std::vector<Result> read_legacy_esd(std::istream &in, std::string &center, LensSummary &sample,
                                    bool companion = false) {
    std::vector<Result> results;
    for (std::string line; std::getline(in, line);) {
        line = bcg::trim(line);
        if (line.empty()) continue;
        if (!companion && line.rfind("# fit_input=", 0) == 0)
            throw std::runtime_error("Use the companion ESD fit input file " + line.substr(12));
        if (line.rfind("# measurement=", 0) == 0)
            throw std::runtime_error("Use the binary measurement snapshot " + line.substr(14));
        if (line.rfind("# center=", 0) == 0 && center.empty()) {
            auto fields = legacy_header(std::move(line));
            read_field(fields, "center", center);
            read_field(fields, "groups", sample.n_lenses);
            read_field(fields, "z_eff", sample.redshift);
            read_field(fields, "mean_Mstar_Msun_over_h", sample.stellar_mass);
            read_field(fields, "n_jack", sample.regions);
            expect_end(fields);
        } else if (line.front() != '#') {
            if (results.size() >= max_radial_bins) throw std::runtime_error("Too many radial bins");
            std::istringstream row(line);
            Result result;
            result.low = table_number(row);
            result.high = table_number(row);
            result.effective = table_number(row);
            if (!(row >> result.count)) throw std::runtime_error("Invalid lens–source pair count");
            result.measured = table_number(row);
            result.sigma = table_number(row);
            result.cross = table_number(row);
            result.cross_sigma = table_number(row);
            expect_end(row);
            results.push_back(std::move(result));
        }
    }
    if (in.bad() || center.empty() || results.empty())
        throw std::runtime_error("Expected a saved ESD measurement with lens sample metadata");
    return results;
}

Covariance read_legacy_covariance(const bcg::fs::path &path, const std::vector<Result> &results,
                                   const std::string &center, int regions) {
    auto in = open_input(path);
    const auto n = results.size();
    Covariance covariance;
    bool have_header = false;
    for (std::string line; std::getline(in, line);) {
        line = bcg::trim(line);
        if (line.empty()) continue;
        if (line.rfind("# center=", 0) == 0 && !have_header) {
            auto fields = legacy_header(std::move(line));
            std::string cov_center;
            int cov_regions = 0;
            read_field(fields, "center", cov_center);
            read_field(fields, "n_jack", cov_regions);
            if (cov_center != center || cov_regions != regions)
                throw std::runtime_error("ESD table and covariance disagree on center or jackknife regions");
            expect(fields, "fit_bins_zero_based");
            while (!(fields >> std::ws).eof()) {
                int bin = 0;
                if (covariance.fit_bins.size() == n || !(fields >> bin))
                    throw std::runtime_error("Invalid fitted bin in covariance header");
                covariance.fit_bins.push_back(bin);
            }
            validate_bins(covariance.fit_bins, static_cast<int>(n));
            have_header = true;
        } else if (line.front() != '#') {
            if (covariance.matrix.size() == n * n) throw std::runtime_error("Too many covariance rows");
            std::istringstream row(line);
            for (std::size_t i = 0; i < n; ++i) covariance.matrix.push_back(table_number(row));
            expect_end(row);
        }
    }
    if (in.bad() || !have_header || covariance.matrix.size() != n * n)
        throw std::runtime_error("Covariance dimensions or fitted bin metadata do not match the ESD table");
    for (const int bin : covariance.fit_bins) {
        const double sigma = results[bin].sigma;
        const double expected = std::sqrt(covariance.matrix[bin * n + bin]);
        const double tolerance = 8 * std::numeric_limits<double>::epsilon() * std::max(sigma, expected);
        if (!(std::abs(sigma - expected) <= tolerance))
            throw std::runtime_error("ESD table uncertainties do not match the covariance diagonal");
    }
    return covariance;
}
} // namespace

FitInput make_fit_input(const std::vector<Result> &results, const Covariance &covariance,
                        const LensSummary &sample, bcg::Center center, const CosmologyParameters &cosmology) {
    FitInput input;
    input.metadata = {center, cosmology, sample, static_cast<int>(results.size())};
    input.bins = covariance.fit_bins;
    for (const int i : input.bins) {
        input.radius.push_back(results[i].effective);
        input.measured.push_back(results[i].measured);
        for (const int j : input.bins)
            input.covariance.push_back(covariance.matrix[i * results.size() + j]);
    }
    return input;
}

namespace {
using bcg::fs::path;

path legacy_snapshot(const path &table) {
    return table.parent_path() / (table.stem().string() + "_fit_input.dat");
}

path input_path(const Config &cfg) {
    if (cfg.stage == Stage::all) return cfg.measurement_output;
    if (!cfg.fit_input_path.empty()) return cfg.fit_input_path;
    for (const auto &candidate : {path(cfg.measurement_output), legacy_snapshot(cfg.legacy_measurement),
                                  legacy_snapshot(cfg.output), path(cfg.legacy_measurement), path(cfg.output)})
        if (bcg::fs::is_regular_file(candidate)) return candidate;
    return cfg.measurement_output;
}

FitInput from_measurement(Measurement &measurement, const path &source, bool diagonal_only) {
    const auto &[center, cosmology, sample, radial_bins] = measurement.metadata;
    const auto covariance = jackknife_covariance(measurement.bins, sample.regions, diagonal_only);
    auto input = make_fit_input(measurement.bins, covariance, sample, center, cosmology);
    input.measurement_id = measurement.id;
    input.measurement_path = source.string();
    return input;
}

bool legacy_fit_snapshot(const path &source) {
    return bcg::archive::has_magic(source, format_name);
}

void hash_file(std::uint64_t &hash, const path &source) {
    std::ifstream in(source, std::ios::binary);
    char buffer[16384];
    while (in.read(buffer, sizeof(buffer)) || in.gcount())
        bcg::archive::hash_bytes(hash, buffer, static_cast<std::size_t>(in.gcount()));
    if (!in.eof() || in.bad()) throw std::runtime_error("Cannot hash measurement input: " + source.string());
}

FitInput load_legacy_fit(const path &source, const Config &cfg) {
    auto in = open_input(source);
    FitInput input;
    auto id = bcg::archive::hash_start;
    hash_file(id, source);
    if (legacy_fit_snapshot(source)) {
        input = read_snapshot(in);
    } else {
        std::string center;
        LensSummary sample;
        const auto bins = read_legacy_esd(in, center, sample);
        const auto cov_file = source.parent_path() / (source.stem().string() + "_cov.dat");
        const auto covariance = read_legacy_covariance(cov_file, bins, center, sample.regions);
        hash_file(id, cov_file);
        std::cout << "Legacy ESD measurement: using [measurement] cosmology from the fit configuration" << std::endl;
        input = make_fit_input(bins, covariance, sample, bcg::parse_center(center), cfg.cosmology);
        input.legacy_cosmology = true;
    }
    input.measurement_id = id;
    input.measurement_path = source.string();
    return input;
}

path companion_table(const path &source, const Config &cfg) {
    const std::string name = source.filename().string(), suffix = "_fit_input.dat";
    if (name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
        const auto prefix = name.substr(0, name.size() - suffix.size());
        const auto filename = path(cfg.output).filename();
        if (filename.stem() == prefix && bcg::fs::is_regular_file(source.parent_path() / filename))
            return source.parent_path() / filename;
        if (bcg::fs::is_regular_file(source.parent_path() / (prefix + ".dat")))
            return source.parent_path() / (prefix + ".dat");
        for (const auto &entry : bcg::fs::directory_iterator(source.parent_path()))
            if (entry.path().stem() == prefix && entry.is_regular_file()) return entry.path();
    }
    throw std::runtime_error("Full ESD table required beside the legacy fit input: " + source.string());
}

void check_companion(const path &table, const FitInput &fit) {
    auto in = open_input(table);
    std::string center;
    LensSummary sample;
    const auto rows = read_legacy_esd(in, center, sample, true);
    const auto &reference = fit.metadata.sample;
    if (center != bcg::center_name(fit.metadata.center) || sample.n_lenses != reference.n_lenses ||
        sample.redshift != reference.redshift || sample.stellar_mass != reference.stellar_mass ||
        sample.regions != reference.regions || rows.size() != static_cast<std::size_t>(fit.metadata.radial_bins))
        throw std::runtime_error("Legacy ESD table has different measurement metadata");
    for (std::size_t i = 0; i < fit.bins.size(); ++i) {
        const auto &row = rows.at(fit.bins[i]);
        const double error = std::sqrt(fit.covariance[i*(fit.bins.size()+1)]);
        if (row.effective != fit.radius[i] || row.measured != fit.measured[i] || row.sigma != error)
            throw std::runtime_error("Legacy ESD table and fit input contain different data points");
    }
}

std::string table_contents(const path &source) {
    auto in = open_input(source);
    std::ostringstream out;
    out << in.rdbuf();
    if (in.bad() || out.fail()) throw std::runtime_error("Cannot read ESD table: " + source.string());
    return out.str();
}

void protect_measurement_file(const path &source, const Config &cfg, bool export_copy = false) {
    const auto protected_path = bcg::fs::weakly_canonical(source);
    for (const auto &destination : {cfg.output, cfg.fit_plot, cfg.posterior_plot}) {
        if (export_copy && destination == cfg.output) continue;
        if (bcg::fs::weakly_canonical(bcg::fs::absolute(destination)) == protected_path)
            throw std::runtime_error("Pull output overlaps a measurement input: " + destination);
    }
}

std::array<std::vector<FdCalibrationPoint>, 2> read_fd_points(std::istream &in, const Measurement &measurement) {
    std::string line;
    if (!std::getline(in, line) || bcg::trim(line) != fd_points_header)
        throw std::runtime_error("Invalid FD point table header");
    const auto source_count = "# all_R_sources=" + std::to_string(measurement.calibration_sources);
    bool matched_count = false;
    std::array<std::vector<FdCalibrationPoint>, 2> components;
    while (std::getline(in, line)) {
        line = bcg::trim(line);
        if (line.empty()) continue;
        if (line.front() == '#') {
            matched_count |= line == source_count;
            continue;
        }
        std::istringstream row(line);
        int component = 0, bin = 0;
        if (!(row >> component >> bin) || component < 1 || component > 2 || bin < 1 || bin > measurement.fd_bins)
            throw std::runtime_error("Invalid FD component or bin number");
        auto &points = components[component - 1];
        const double low = table_number(row), high = table_number(row);
        FdCalibrationPoint point;
        point.field = table_number(row);
        point.value = table_number(row);
        point.sigma = table_number(row);
        if (!(row >> point.count >> point.status) || point.count < 0 || point.count > measurement.calibration_sources ||
            point.status < 0 || point.status > 2 || static_cast<std::size_t>(bin) != points.size() + 1 ||
            !(low >= -measurement.fd_limit && high <= measurement.fd_limit && low < high))
            throw std::runtime_error("Invalid FD point row");
        if (point.status == 0 && !(point.count > 0 && std::isfinite(point.field) &&
                std::isfinite(point.value) && std::isfinite(point.sigma) && point.sigma > 0))
            throw std::runtime_error("Invalid FD shear estimate or uncertainty");
        expect_end(row);
        points.push_back(point);
    }
    const auto count = static_cast<std::size_t>(measurement.fd_bins);
    if (!in.eof() || in.bad() || !matched_count || components[0].size() != count || components[1].size() != count)
        throw std::runtime_error("FD point table does not match the measurement source or bin counts");
    return components;
}

FdPlot load_fd_plot(const path &snapshot, const Measurement &measurement, const Config &cfg) {
    if (!measurement.fd_bins) return {};
    const std::string suffix = "_measurement";
    auto prefix = snapshot.stem().string();
    if (prefix.size() > suffix.size() && prefix.compare(prefix.size() - suffix.size(), suffix.size(), suffix) == 0)
        prefix.resize(prefix.size() - suffix.size());
    const auto source = snapshot.parent_path() / (prefix + "_fd_points.dat");
    if (!bcg::fs::is_regular_file(source)) {
        std::cout << "FD plot unavailable: missing point table " << source.string() << std::endl;
        return {};
    }
    protect_measurement_file(source, cfg);

    FdPlot result;
    result.calibration = measurement.calibration;
    result.limit = measurement.fd_limit;
    const auto &cal = result.calibration;
    try {
        auto in = open_input(source);
        const auto points = read_fd_points(in, measurement);
        const std::array<double, 2> m{cal.m1, cal.m2}, c{cal.c1, cal.c2};
        for (std::size_t i = 0; i < points.size(); ++i) {
            // Check the companion table against the calibration used for this measurement.
            const auto fit = fit_selfcal_line(points[i]);
            if (std::abs(fit.m - m[i]) > 1e-10 * std::max(1.0, std::abs(m[i])) ||
                std::abs(fit.c - c[i]) > 1e-10 * std::max(result.limit, std::abs(c[i])))
                throw std::runtime_error("FD points do not reproduce the saved calibration");
            for (const auto &p : points[i])
                if (p.status == 0) result.points[i].insert(result.points[i].end(), {p.field, p.value, p.sigma});
        }
    } catch (const std::runtime_error &error) {
        std::cerr << "FD plot unavailable for " << source.string() << ": " << error.what() << std::endl;
        return {};
    }

    result.output = cfg.legacy_fit_base + "_fd.png";
    const auto destination = bcg::fs::weakly_canonical(bcg::fs::absolute(result.output));
    for (const auto &reserved : {cfg.output, cfg.fit_plot, cfg.posterior_plot, cfg.fit_archive,
                                 cfg.measurement_output, snapshot.string(), source.string()}) {
        if (bcg::fs::weakly_canonical(bcg::fs::absolute(reserved)) == destination)
            throw std::runtime_error("FD plot output must have a distinct path: " + result.output);
    }
    return result;
}
} // namespace

FitInput load_fit_input(const Config &cfg) {
    const auto source = bcg::fs::weakly_canonical(bcg::fs::absolute(input_path(cfg)));
    std::cout << "Reading ESD fit input: " << source.string() << std::endl;
    if (is_measurement_snapshot(source) || source.extension() == ".bin") {
        auto measurement = read_measurement(source);
        return from_measurement(measurement, source, false);
    }
    try {
        return load_legacy_fit(source, cfg);
    } catch (const std::exception &error) {
        throw std::runtime_error("Invalid ESD fit input " + source.string() + ": " + error.what());
    }
}

PullInput load_pull_input(const Config &cfg, const std::string &saved_measurement, bool legacy_fit) {
    const auto requested = (cfg.stage == Stage::pull && !cfg.fit_input_path.empty()) || saved_measurement.empty()
        ? input_path(cfg) : path(saved_measurement);
    const auto source = bcg::fs::weakly_canonical(bcg::fs::absolute(requested));
    std::cout << "Reading ESD measurement for pull: " << source.string() << std::endl;
    if (is_measurement_snapshot(source) || source.extension() == ".bin") {
        protect_measurement_file(source, cfg);
        auto measurement = read_measurement(source);
        auto fit = from_measurement(measurement, source, !legacy_fit);
        validate_fit_input(fit, !legacy_fit);
        if (cfg.stage == Stage::pull) print_measurement_summary(measurement);
        return {std::move(fit), format_esd_table(measurement, source.filename().string()),
                load_fd_plot(source, measurement, cfg)};
    }
    auto fit = load_legacy_fit(source, cfg);
    validate_fit_input(fit);
    const bool companion = legacy_fit_snapshot(source);
    const auto table = companion ? companion_table(source, cfg) : source;
    protect_measurement_file(source, cfg, !companion);
    if (companion) {
        protect_measurement_file(table, cfg, true);
        check_companion(table, fit);
    } else {
        protect_measurement_file(source.parent_path() / (source.stem().string() + "_cov.dat"), cfg);
    }
    return {std::move(fit), table_contents(table), {}};
}

void validate_fit_input(const FitInput &input, bool diagonal_only) {
    const auto &[center, cosmology, sample, radial_bins] = input.metadata;
    validate_cosmology(cosmology);
    if (sample.regions < 2 || sample.n_lenses < sample.regions ||
        !(std::isfinite(sample.redshift) && sample.redshift > 0) || !std::isfinite(sample.stellar_mass) ||
        (center == bcg::Center::bcg ? sample.stellar_mass <= 0 : sample.stellar_mass != 0))
        throw std::runtime_error("Invalid lens sample in ESD fit input");
    const auto n = input.bins.size();
    if (radial_bins < 1 || radial_bins > max_radial_bins || n > static_cast<std::size_t>(radial_bins) ||
        input.radius.size() != n || input.measured.size() != n || input.covariance.size() != n * n)
        throw std::runtime_error("Invalid dimensions in ESD fit input");
    validate_bins(input.bins, radial_bins);
    for (std::size_t i = 0; i < n; ++i) {
        if (!(input.radius[i] > 0 && std::isfinite(input.radius[i]) && std::isfinite(input.measured[i])))
            throw std::runtime_error("Invalid fitted radius or ESD in measurement");
        for (std::size_t j = diagonal_only ? i : 0; j <= i; ++j) {
            const double value = input.covariance[i * n + j];
            if (!std::isfinite(value) || value != input.covariance[j * n + i])
                throw std::runtime_error("Fitted ESD covariance must be finite and symmetric");
        }
    }
}
} // namespace bcg_esd
