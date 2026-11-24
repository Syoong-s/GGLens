#include "bcg_esd.hpp"
#include <iostream>
#include <iterator>
#include <variant>

namespace bcg_esd {
Stage parse_stage(const std::string &name) {
    if (name == "all") return Stage::all;
    if (name == "measure") return Stage::measure;
    if (name == "fit") return Stage::fit;
    if (name == "pull") return Stage::pull;
    throw std::invalid_argument("stage must be all, measure, fit or pull");
}

void validate_cosmology(const CosmologyParameters &c) {
    if (!(c.omega_m > 0 && c.omega_m < 1 && std::isfinite(c.h) && c.h > 0 &&
          std::isfinite(c.sigma8) && c.sigma8 > 0 && std::isfinite(c.spectral_index) && c.spectral_index > 0 &&
          c.omega_b_h2 > 0 && c.neutrino_mass >= 0 &&
          c.omega_m*c.h*c.h > c.omega_b_h2 + c.neutrino_mass/93.14))
        throw std::runtime_error("Invalid cosmology parameters");
}

namespace {
// ==========================================
// Function: parse_ccd_cut
// Method: Read a comma-separated list of CCD IDs or the common disabled value.
// ==========================================
std::vector<int> parse_ccd_cut(const std::string &input) {
    const std::string value = bcg::trim(input);
    if (value.empty()) return {};
    if (value.front() == '-' && value.find(',') == std::string::npos) {
        std::size_t end = 0;
        const double threshold = std::stod(value, &end);
        if (end == value.size() && std::isfinite(threshold) && !bcg::cut_enabled(threshold)) return {};
    }

    if (value.back() == ',') throw std::invalid_argument("expected a CCD ID after the last comma");
    std::vector<int> ids;
    std::istringstream entries(value);
    for (std::string token; std::getline(entries, token, ',');) {
        token = bcg::trim(token);
        std::size_t end = 0;
        const int id = std::stoi(token, &end);
        if (end != token.size() || id < 0)
            throw std::invalid_argument("expected nonnegative integer CCD IDs separated by commas");
        ids.push_back(id);
    }
    return ids;
}

// ==========================================
// Function: apply_ini
// Method: Overlay typed INI keys on compiled defaults and report invalid keys at their line.
// ==========================================
void apply_ini(Config &c, const std::string &text) {
    using Setting = std::variant<std::string *, double *, int *, bool *, Stage *, bcg::Center *, std::vector<int> *,
                                 bcg::CatalogDelimiter *, bcg::MassFormat *>;
    auto &[group, galaxy, member] = c.lens_columns;
    const std::unordered_map<std::string, Setting> settings{
        {"run.stage", &c.stage},
        {"run.group-list", &c.group_list},
        {"run.shear-list", &c.shear_list},
        {"run.output", &c.output},
        {"run.output-dir", &c.output_dir},
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
        {"measurement.omega-m", &c.cosmology.omega_m},
        {"measurement.h", &c.cosmology.h},
        {"measurement.sigma8", &c.cosmology.sigma8},
        {"measurement.omega-b-h2", &c.cosmology.omega_b_h2},
        {"measurement.spectral-index", &c.cosmology.spectral_index},
        {"measurement.neutrino-mass", &c.cosmology.neutrino_mass},
        {"measurement.n-r-bins", &c.n_r_bins},
        {"measurement.n-jack", &c.n_jack},
        {"measurement.pdf-bins", &c.pdf_bins},
        {"group_columns.id", &group.id},
        {"group_columns.redshift", &group.redshift},
        {"group_columns.log-halo-mass", &group.halo_mass},
        {"group_columns.delimiter", &group.text.delimiter},
        {"group_columns.header", &group.text.header},
        {"group_columns.mass-format", &group.mass_format},
        {"galaxy_columns.id", &galaxy.id},
        {"galaxy_columns.ra", &galaxy.ra},
        {"galaxy_columns.dec", &galaxy.dec},
        {"galaxy_columns.luminosity", &galaxy.luminosity},
        {"galaxy_columns.log-stellar-mass", &galaxy.stellar_mass},
        {"galaxy_columns.photometry-select", &galaxy.photometry_select},
        {"galaxy_columns.delimiter", &galaxy.text.delimiter},
        {"galaxy_columns.header", &galaxy.text.header},
        {"galaxy_columns.mass-format", &galaxy.mass_format},
        {"member_columns.galaxy-id", &member.galaxy_id},
        {"member_columns.group-id", &member.group_id},
        {"member_columns.brightest-flag", &member.brightest},
        {"member_columns.brightest-value", &member.brightest_value},
        {"member_columns.delimiter", &member.text.delimiter},
        {"member_columns.header", &member.text.header},
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
        {"external_columns.ccd-num", &c.external_columns.ccd_num},
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
        {"cuts.ccd-cut", &c.cuts.excluded_ccds},
        {"fit.input", &c.fit_input_path},
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
        "run", "measurement", "group_columns", "galaxy_columns", "member_columns",
        "calibration", "external_columns", "cuts", "fit", "plots"
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
                } else if constexpr (std::is_same_v<T, std::vector<int>>) {
                    *dest = parse_ccd_cut(value);
                } else if constexpr (std::is_same_v<T, Stage>) {
                    *dest = parse_stage(value);
                } else if constexpr (std::is_same_v<T, bcg::Center>) {
                    *dest = bcg::parse_center(value);
                } else if constexpr (std::is_same_v<T, bcg::CatalogDelimiter>) {
                    if (value == "whitespace") *dest = bcg::CatalogDelimiter::whitespace;
                    else if (value == "comma") *dest = bcg::CatalogDelimiter::comma;
                    else if (value == "tab") *dest = bcg::CatalogDelimiter::tab;
                    else throw std::invalid_argument("expected whitespace, comma or tab");
                } else if constexpr (std::is_same_v<T, bcg::MassFormat>) {
                    if (value == "log10") *dest = bcg::MassFormat::log10;
                    else if (value == "linear") *dest = bcg::MassFormat::linear;
                    else throw std::invalid_argument("expected log10 or linear");
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
// Function: validate_measurement_config
// Method: Check the settings used by catalog ingestion and measurement.
// ==========================================
void validate_measurement_config(Config &c) {
    if (!c.self_calibrate) c.calibration.validate();
    if (c.group_list.empty() || c.shear_list.empty())
        throw std::runtime_error("run.group-list and run.shear-list must be set for measurement");
    validate_cosmology(c.cosmology);
    if (!(c.zmin >= 0 && c.zmax > c.zmin && c.mmax > c.mmin && c.rmin > 0 && c.rmax > c.rmin &&
          c.dzmin >= 0 && c.n_r_bins >= 1 && c.n_jack >= 2 && c.pdf_bins > 0))
        throw std::runtime_error("Invalid measurement parameter range");
    if (c.stage == Stage::all && (c.n_r_bins < 6 || c.n_jack < 7))
        throw std::runtime_error("A combined run needs at least six radial bins and seven requested jackknife regions");
    auto &excluded_ccds = c.cuts.excluded_ccds;
    std::sort(excluded_ccds.begin(), excluded_ccds.end());
    excluded_ccds.erase(std::unique(excluded_ccds.begin(), excluded_ccds.end()), excluded_ccds.end());
    if (!excluded_ccds.empty() && excluded_ccds.front() < 0)
        throw std::runtime_error("cuts.ccd-cut requires nonnegative CCD IDs");
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
        static_cast<std::uint64_t>(c.n_r_bins) * 16384 > budget / 8 || c.n_r_bins > max_radial_bins)
        throw std::runtime_error("Bin/replica/duplicate configuration exceeds memory or open-file budget");
    if (c.self_calibrate) {
        if (!(c.fd_limit() > 0 && c.selfcal_fd_bins >= 2 && c.selfcal_fd_bins <= 100 &&
              c.selfcal_pdf_bins >= 2 && c.selfcal_pdf_bins <= 200 &&
              c.selfcal_min_count > c.selfcal_pdf_bins && c.selfcal_search_width > 0))
            throw std::runtime_error("Invalid self-calibration FD settings");
    }
}

void validate_fit_config(const Config &c) {
    if (c.fit_walkers < 10 || c.fit_walkers % 2 || c.fit_steps <= c.fit_burn ||
        c.fit_burn < 0 || c.fit_ensembles < 1 || c.fit_seed < 0)
        throw std::runtime_error("Invalid MCMC settings: use >=10 even walkers and steps > burn-in");
    if (c.fit_profile_draws < 1)
        throw std::runtime_error("fit.profile-draws must be positive");
    const auto chain_values = static_cast<long double>(c.fit_walkers) *
                              (c.fit_steps - c.fit_burn) * c.fit_ensembles * 6;
    if (chain_values > std::numeric_limits<int>::max())
        throw std::runtime_error("MCMC chain exceeds MPI count limits");
}

void resolve_outputs(Config &c, int rank) {
    using bcg::fs::path;
    const path output(c.output), directory(c.output_dir);
    if (output.filename().empty() || output.filename() == "." || output.filename() == "..")
        throw std::runtime_error("run.output must name an ESD file");
    if (!directory.empty() && output.has_parent_path())
        throw std::runtime_error("run.output must be a filename when run.output-dir is set");

    const path result_root = directory.empty() ? output.parent_path() : directory;
    const auto filename = output.filename();
    const auto stem = output.stem().string();
    if (filename == "measurement" || filename == "fit")
        throw std::runtime_error("run.output conflicts with a stage directory: " + filename.string());
    c.output = (result_root / filename).string();
    c.measurement_output = (result_root / "measurement" / (stem + "_measurement.bin")).string();
    c.legacy_measurement = (result_root / "measurement" / filename).string();
    c.fit_archive = (result_root / "fit" / stem).string();
    c.legacy_fit_base = (result_root / stem).string();

    const bool measures = c.stage == Stage::all || c.stage == Stage::measure;
    const bool renders = c.stage == Stage::all || c.stage == Stage::pull;
    if (measures) {
        if (!c.cache.empty() && c.center != bcg::Center::bcg) {
            const path cache(c.cache);
            c.cache = (cache.parent_path() / (cache.stem().string() + "_" +
                       bcg::center_name(c.center) + cache.extension().string())).string();
        }
        if (c.scratch_dir.empty())
            c.scratch_dir = bcg::fs::absolute(result_root / ".bcg_scratch").string();
    }
    if (renders) {
        c.fit_plot = c.fit_plot.empty() ? c.legacy_fit_base + "_fit.png" : (directory / c.fit_plot).string();
        c.posterior_plot = c.legacy_fit_base + "_posterior.png";
        if (path(c.fit_plot).extension() != ".png")
            throw std::runtime_error("Plot output must end in .png: " + c.fit_plot);
    }

    if (rank == 0) {
        const auto canonical = [](const std::string &file) {
            return bcg::fs::weakly_canonical(bcg::fs::absolute(file)).string();
        };
        std::unordered_set<std::string> paths{canonical(c.measurement_output), canonical(c.fit_archive)};
        const auto require_distinct = [&](const std::string &file) {
            if (!file.empty() && !paths.insert(canonical(file)).second)
                throw std::runtime_error("Each output must have a distinct path: " + file);
        };
        if (renders) {
            require_distinct(c.output);
            require_distinct(c.fit_plot);
            require_distinct(c.posterior_plot);
        }
        if (measures) {
            if (!renders) require_distinct(c.output);
            require_distinct(c.cache);
        }
    }
}
} // namespace

// ==========================================
// Function: load_config
// Method: Read the supplied INI or the default file on rank zero and broadcast its text.
// ==========================================
Config load_config(const char *ini_path, MPI_Comm comm, const char *stage_override) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    std::string text;
    if (rank == 0) {
        const bool explicit_path = ini_path != nullptr;
        const std::string path = explicit_path ? ini_path : defaults::ini_name;
        if (explicit_path || bcg::fs::exists(path)) {
            if (!bcg::fs::is_regular_file(path))
                throw std::runtime_error("INI is not a regular file: '" + path + "'");
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
    if (stage_override) cfg.stage = parse_stage(stage_override);
    if (cfg.stage == Stage::measure || cfg.stage == Stage::all) validate_measurement_config(cfg);
    if (cfg.stage == Stage::fit || cfg.stage == Stage::all) validate_fit_config(cfg);
    resolve_outputs(cfg, rank);
    return cfg;
}
} // namespace bcg_esd
