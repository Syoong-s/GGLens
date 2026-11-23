#include "bcg_esd.hpp"
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

// Legacy tables use nan for empty bins. Formatted double extraction does not read it.
double table_number(std::istream &in) {
    std::string token;
    if (!(in >> token)) throw std::runtime_error("Incomplete measurement or covariance row");
    char *end = nullptr;
    const double value = std::strtod(token.c_str(), &end);
    if (end != token.c_str() + token.size()) throw std::runtime_error("Invalid table value: " + token);
    return value;
}

std::vector<Result> read_legacy_esd(std::istream &in, std::string &center, LensSummary &sample) {
    std::vector<Result> results;
    for (std::string line; std::getline(in, line);) {
        line = bcg::trim(line);
        if (line.empty()) continue;
        if (line.rfind("# fit_input=", 0) == 0)
            throw std::runtime_error("Use the companion ESD fit input file " + line.substr(12));
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

void write_fit_input(std::ostream &out, const FitInput &input) {
    const auto &[center, cosmology, sample, radial_bins] = input.metadata;
    const auto count = input.bins.size();
    out << format_name << ' ' << format_version << '\n'
        << "center " << bcg::center_name(center) << '\n'
        << "sample " << sample.n_lenses << ' ' << sample.regions << ' ' << sample.redshift << ' ' << sample.stellar_mass << '\n'
        << "cosmology " << cosmology.omega_m << ' ' << cosmology.h << ' ' << cosmology.sigma8 << ' '
        << cosmology.omega_b_h2 << ' ' << cosmology.spectral_index << ' ' << cosmology.neutrino_mass << '\n'
        << "radial_bins " << radial_bins << '\n' << "fit_bins " << count << '\n'
        << "index_radius_esd\n";
    for (std::size_t i = 0; i < count; ++i)
        out << input.bins[i] << ' ' << input.radius[i] << ' ' << input.measured[i] << '\n';
    out << "covariance\n";
    for (std::size_t row = 0; row < count; ++row)
        for (std::size_t col = 0; col < count; ++col)
            out << input.covariance[row * count + col] << (col + 1 == count ? '\n' : ' ');
}

FitInput load_fit_input(const Config &cfg) {
    bcg::fs::path path(cfg.fit_input_path.empty() ? cfg.fit_input_output : cfg.fit_input_path);
    if (cfg.fit_input_path.empty() && !bcg::fs::exists(path) && bcg::fs::exists(cfg.output))
        path = cfg.output;
    auto in = open_input(path);
    std::cout << "Reading ESD fit input: " << path.string() << std::endl;
    try {
        std::string first;
        in >> first;
        in.seekg(0);
        if (first == format_name) return read_snapshot(in);

        std::string center;
        LensSummary sample;
        const auto results = read_legacy_esd(in, center, sample);
        const auto cov_path = path.parent_path() / (path.stem().string() + "_cov.dat");
        const auto covariance = read_legacy_covariance(cov_path, results, center, sample.regions);
        std::cout << "Legacy ESD measurement: using [measurement] cosmology from the fit configuration" << std::endl;
        return make_fit_input(results, covariance, sample, bcg::parse_center(center), cfg.cosmology);
    } catch (const std::exception &error) {
        throw std::runtime_error("Invalid ESD fit input " + path.string() + ": " + error.what());
    }
}

void validate_fit_input(const FitInput &input) {
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
        for (std::size_t j = 0; j <= i; ++j) {
            const double value = input.covariance[i * n + j];
            if (!std::isfinite(value) || value != input.covariance[j * n + i])
                throw std::runtime_error("Fitted ESD covariance must be finite and symmetric");
        }
    }
}
} // namespace bcg_esd
