#include "bcg_esd.hpp"
#include <iostream>
#include <iterator>
#include <variant>

namespace bcg_esd {
namespace {
// ==========================================
// Function: apply_ini
// Method: Overlay typed INI keys on compiled defaults and report invalid keys at their line.
// ==========================================
void apply_ini(Config &c, const std::string &text) {
    using Setting = std::variant<std::string *, double *, int *, bool *, bcg::Center *>;
    const std::unordered_map<std::string, Setting> settings{
        {"run.group-list", &c.group_list},
        {"run.shear-list", &c.shear_list},
        {"run.output", &c.output},
        {"run.cache", &c.cache},
        {"run.scratch-dir", &c.scratch_dir},
        {"run.memory-mb", &c.memory_mib},
        {"measurement.center", &c.center},
        {"measurement.z-min", &c.zmin},
        {"measurement.z-max", &c.zmax},
        {"measurement.logm-min", &c.mmin},
        {"measurement.logm-max", &c.mmax},
        {"measurement.r-min", &c.rmin},
        {"measurement.r-max", &c.rmax},
        {"measurement.dz-min", &c.dzmin},
        {"measurement.omega-m", &c.omega_m},
        {"measurement.h", &c.h},
        {"measurement.sigma8", &c.sigma8},
        {"measurement.omega-b-h2", &c.omega_b_h2},
        {"measurement.spectral-index", &c.spectral_index},
        {"measurement.neutrino-mass", &c.neutrino_mass},
        {"measurement.n-r-bins", &c.n_r_bins},
        {"measurement.n-jack", &c.n_jack},
        {"measurement.pdf-bins", &c.pdf_bins},
        {"calibration.fd-corrected", &c.fd_corrected},
        {"calibration.self-calibrate", &c.self_calibrate},
        {"calibration.m1", &c.calibration.m1},
        {"calibration.m2", &c.calibration.m2},
        {"calibration.c1", &c.calibration.c1},
        {"calibration.c2", &c.calibration.c2},
        {"calibration.gf-lim", &c.cuts.gf_lim},
        {"calibration.selfcal-fd-bins", &c.selfcal_fd_bins},
        {"calibration.selfcal-pdf-bins", &c.selfcal_pdf_bins},
        {"calibration.selfcal-min-count", &c.selfcal_min_count},
        {"calibration.selfcal-search-width", &c.selfcal_search_width},
        {"external_columns.ra", &c.external_columns.ra},
        {"external_columns.dec", &c.external_columns.dec},
        {"external_columns.zp", &c.external_columns.zp},
        {"external_columns.mag", &c.external_columns.mag},
        {"external_columns.ext-mash", &c.external_columns.ext_mash},
        {"external_columns.expo-num", &c.external_columns.expo_num},
        {"external_columns.expo-chi2", &c.external_columns.expo_chi2},
        {"cuts.snrfcut", &c.cuts.snrfcut},
        {"cuts.snrlow", &c.cuts.snrlow},
        {"cuts.snrhigh", &c.cuts.snrhigh},
        {"cuts.starcut", &c.cuts.starcut},
        {"cuts.chi2-thresh", &c.cuts.chi2_thresh},
        {"cuts.flagcut", &c.cuts.flagcut},
        {"cuts.zplow", &c.cuts.zplow},
        {"cuts.zphigh", &c.cuts.zphigh},
        {"cuts.imaxcut", &c.cuts.imaxcut},
        {"cuts.jmaxcut", &c.cuts.jmaxcut},
        {"cuts.mag-low", &c.cuts.mag_low},
        {"cuts.mag-high", &c.cuts.mag_high},
        {"cuts.psf-chi2-mltp", &c.cuts.psf_chi2_mltp},
        {"cuts.r-half-thresh", &c.cuts.r_half_thresh},
        {"cuts.ext-cut", &c.cuts.ext_cut},
        {"cuts.max-dup", &c.cuts.max_dup},
        {"cuts.star-bar-mltp", &c.cuts.star_bar_mltp},
        {"cuts.expo-bar-mltp", &c.cuts.expo_bar_mltp},
        {"fit.walkers", &c.fit_walkers},
        {"fit.steps", &c.fit_steps},
        {"fit.burn-in", &c.fit_burn},
        {"fit.ensembles", &c.fit_ensembles},
        {"fit.seed", &c.fit_seed},
        {"fit.profile-draws", &c.fit_profile_draws},
        {"fit.python-path", &c.fit_python_path},
        {"plots.esd-plot", &c.fit_plot},
        {"plots.title", &c.plot_title},
    };
    const std::unordered_set<std::string> sections{
        "run", "measurement", "calibration", "external_columns", "cuts", "fit", "plots"
    };
    std::unordered_set<std::string> seen;
    std::istringstream input(text);
    std::string section, line;
    int lineno = 0;
    while (std::getline(input, line)) {
        ++lineno;
        line = bcg::trim(line);
        if (line.empty() || line.front() == '#' || line.front() == ';') continue;
        const auto error = [lineno](const std::string &message) {
            return std::runtime_error("INI line " + std::to_string(lineno) + ": " + message);
        };
        if (line.front() == '[' && line.back() == ']') {
            section = bcg::trim(line.substr(1, line.size() - 2));
            if (!sections.count(section)) throw error("Unknown section " + section);
            continue;
        }
        const auto equal = line.find('=');
        if (equal == std::string::npos) throw error("Expected key=value");
        const std::string key = section + "." + bcg::trim(line.substr(0, equal));
        const auto found = settings.find(key);
        if (found == settings.end()) throw error("Unknown setting " + key);
        if (!seen.insert(key).second) throw error("Repeated setting " + key);
        std::string value = bcg::trim(line.substr(equal + 1));
        if (value.size() >= 2 && ((value.front() == '"' && value.back() == '"') ||
                                  (value.front() == '\'' && value.back() == '\'')))
            value = value.substr(1, value.size() - 2);
        try {
            std::visit([&](auto dest) {
                using T = std::remove_pointer_t<decltype(dest)>;
                std::size_t end = 0;
                if constexpr (std::is_same_v<T, std::string>) {
                    *dest = value;
                } else if constexpr (std::is_same_v<T, bcg::Center>) {
                    *dest = bcg::parse_center(value);
                } else if constexpr (std::is_same_v<T, bool>) {
                    if (value == "true" || value == "1") *dest = true;
                    else if (value == "false" || value == "0") *dest = false;
                    else throw std::invalid_argument("expected true/false or 1/0");
                } else {
                    if constexpr (std::is_same_v<T, double>) *dest = std::stod(value, &end);
                    else *dest = std::stoi(value, &end);
                    if (end != value.size() || !std::isfinite(*dest))
                        throw std::invalid_argument("expected a finite number");
                }
            }, found->second);
        } catch (const std::logic_error &e) {
            throw error("Invalid value for " + key + ": " + value + " (" + e.what() + ")");
        }
    }
}

// ==========================================
// Function: validate_config
// Method: Validate ranges and resolve output defaults after the INI overlay.
// ==========================================
void validate_config(Config &c) {
    if (!c.self_calibrate) c.calibration.validate();
    if (c.group_list.empty() || c.shear_list.empty() || c.output.empty())
        throw std::runtime_error("run.group-list, run.shear-list and run.output must be set");
    if (!(c.zmin >= 0 && c.zmax > c.zmin && c.mmax > c.mmin && c.rmin > 0 && c.rmax > c.rmin &&
          c.dzmin >= 0 && c.omega_m > 0 && c.omega_m < 1 && c.h > 0 && c.sigma8 > 0 &&
          c.omega_b_h2 > 0 && c.neutrino_mass >= 0 && c.spectral_index > 0 &&
          c.omega_m*c.h*c.h > c.omega_b_h2 + c.neutrino_mass/93.14 &&
          c.n_r_bins >= 6 && c.n_jack > 6 && c.pdf_bins > 0))
        throw std::runtime_error("Invalid measurement parameter range");
    const auto &cuts = c.cuts;
    using bcg::cut_enabled;
    const bool limit_duplicates = cut_enabled(cuts.max_dup);
    if (limit_duplicates && (!(cuts.max_dup >= 1 && cuts.max_dup <= std::numeric_limits<int>::max()) ||
                             cuts.max_dup != std::trunc(cuts.max_dup)))
        throw std::runtime_error("cuts.max-dup must be a positive integer or below -999");
    if (cut_enabled(cuts.gf_lim) && cuts.gf_lim < 0)
        throw std::runtime_error("calibration.gf-lim must be nonnegative or below -999");
    if ((cut_enabled(cuts.mag_low) && cut_enabled(cuts.mag_high) && cuts.mag_high < cuts.mag_low) ||
        (cut_enabled(cuts.zplow) && cut_enabled(cuts.zphigh) && cuts.zphigh <= cuts.zplow) ||
        (cut_enabled(cuts.snrlow) && cut_enabled(cuts.snrhigh) && cuts.snrhigh < cuts.snrlow))
        throw std::runtime_error("Invalid shear-cut parameter range");
    if (static_cast<std::uint64_t>(c.pdf_bins + 1LL) * 400 > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        c.n_jack == std::numeric_limits<int>::max())
        throw std::runtime_error("PDF/replica dimensions exceed MPI count limits");
    if (c.memory_mib < 64) throw std::runtime_error("run.memory-mb must be at least 64 MiB");
    const auto budget = c.memory_bytes();
    if ((limit_duplicates && static_cast<std::uint64_t>(cuts.max_dup) * sizeof(bcg::ShearRow) > budget / 16) ||
        static_cast<std::uint64_t>(c.pdf_bins + 1LL) * 3200 > budget / 16 ||
        static_cast<std::uint64_t>(c.n_jack + 1LL) * 2049 * sizeof(std::int64_t) > budget / 4 ||
        static_cast<std::uint64_t>(c.n_jack + 1LL) * (c.pdf_bins + 1LL) * 4 * 61 * sizeof(std::int64_t)
            > budget / 4 ||
        static_cast<std::uint64_t>(c.n_r_bins) * 16384 > budget / 8 || c.n_r_bins > 256)
        throw std::runtime_error("Bin/replica/duplicate configuration exceeds memory or open-file budget");
    if (c.self_calibrate) {
        if (!(c.fd_limit() > 0 && c.selfcal_fd_bins >= 2 && c.selfcal_fd_bins <= 100 &&
              c.selfcal_pdf_bins >= 2 && c.selfcal_pdf_bins <= 200 &&
              c.selfcal_min_count > c.selfcal_pdf_bins && c.selfcal_search_width > 0))
            throw std::runtime_error("Invalid self-calibration FD settings");
    }
    if (c.fit_walkers < 10 || c.fit_walkers % 2 || c.fit_steps <= c.fit_burn ||
        c.fit_burn < 0 || c.fit_ensembles < 1 || c.fit_seed < 0)
        throw std::runtime_error("Invalid MCMC settings: use >=10 even walkers and steps > burn-in");
    if (c.fit_profile_draws < 1)
        throw std::runtime_error("fit.profile-draws must be positive");
    const auto chain_values = static_cast<long double>(c.fit_walkers) *
                              (c.fit_steps - c.fit_burn) * c.fit_ensembles * 6;
    if (chain_values > std::numeric_limits<int>::max())
        throw std::runtime_error("MCMC chain exceeds MPI count limits");
    const auto path = bcg::fs::path(c.output);
    const auto base = (path.parent_path() / path.stem()).string();
    c.covariance_output = base + "_cov.dat";
    c.replica_output = base + "_jackknife.dat";
    c.summary_output = base + "_fit.json";
    c.chain_output = base + "_chain.npz";
    c.profile_output = base + "_profiles.dat";
    c.profile_draws_output = base + "_profile_draws.npz";
    if (!c.cache.empty() && c.center != bcg::Center::bcg) {
        const bcg::fs::path cache(c.cache);
        c.cache = (cache.parent_path() / (cache.stem().string() + "_" +
                   bcg::center_name(c.center) + cache.extension().string())).string();
    }
    if (c.plot_title.empty()) {
        std::string label = bcg::center_name(c.center);
        for (char &letter : label) letter = static_cast<char>(std::toupper(static_cast<unsigned char>(letter)));
        c.plot_title = label + "-centered group lensing";
    }
    if (c.fit_plot.empty()) c.fit_plot = base + "_fit.png";
    if (bcg::fs::path(c.fit_plot).extension() != ".png")
        throw std::runtime_error("Plot output must end in .png: " + c.fit_plot);
    std::unordered_set<std::string> outputs;
    for (const auto &output : {c.output, c.cache, c.fit_plot, c.covariance_output,
                               c.replica_output, c.summary_output, c.chain_output,
                               c.profile_output, c.profile_draws_output}) {
        if (!output.empty() && !outputs.insert(bcg::fs::weakly_canonical(
                bcg::fs::absolute(output)).string()).second)
            throw std::runtime_error("Each output must have a distinct path: " + output);
    }
    if (c.scratch_dir.empty()) c.scratch_dir = (bcg::fs::absolute(c.output).parent_path() / ".bcg_scratch").string();
}
} // namespace

// ==========================================
// Function: load_config
// Method: Read the optional override on rank zero and apply the same text on every rank.
// ==========================================
Config load_config(MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    std::string text;
    if (rank == 0) {
        const char *selected = std::getenv(defaults::ini_environment);
        const bool explicit_path = selected && *selected;
        const std::string path = explicit_path ? selected : defaults::ini_name;
        if (explicit_path || bcg::fs::exists(path)) {
            if (!bcg::fs::is_regular_file(path))
                throw std::runtime_error("INI is not a regular file: " + path);
            std::ifstream in(path);
            if (!in) throw std::runtime_error("Cannot open INI: " + path);
            std::ostringstream contents;
            contents << in.rdbuf();
            if (in.bad() || contents.bad()) throw std::runtime_error("Cannot read INI: " + path);
            text = std::move(contents).str();
            std::cout << "INI overrides: " << path << std::endl;
        }
        if (text.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("INI exceeds MPI broadcast limit");
    }
    int size = static_cast<int>(text.size());
    MPI_Bcast(&size, 1, MPI_INT, 0, comm);
    if (rank != 0) text.resize(static_cast<std::size_t>(size));
    MPI_Bcast(text.data(), size, MPI_CHAR, 0, comm);
    Config cfg;
    apply_ini(cfg, text);
    validate_config(cfg);
    return cfg;
}
} // namespace bcg_esd
