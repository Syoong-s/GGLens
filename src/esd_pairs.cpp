#include "bcg_esd.hpp"
#include "bcg_star_bar.hpp"
#include "bcg_physics.hpp"
#include "selfcalibration.hpp"
#include <iostream>

namespace bcg_esd {
// ==========================================
// Function: radial_edges
// Method: Retain the logarithmic h^-1 Mpc edges and half-open radial cuts.
// ==========================================
std::vector<double> radial_edges(const Config &cfg) {
    std::vector<double> edges(static_cast<std::size_t>(cfg.n_r_bins) + 1);
    for (int i = 0; i <= cfg.n_r_bins; ++i)
        edges[i] = std::pow(10.0, std::log10(cfg.rmin) +
            (std::log10(cfg.rmax) - std::log10(cfg.rmin)) * i / cfg.n_r_bins);
    edges.front() = cfg.rmin; edges.back() = cfg.rmax; return edges;
}

// ==========================================
// Function: build_pair_spools
// Method: Restore catalog FD when enabled, fit raw all-R sources, and spool calibrated ESD pairs.
// ==========================================
std::vector<std::unique_ptr<bcg::PairStore>> build_pair_spools(
    const std::vector<bcg::Lens> &lenses, const std::vector<int> &labels,
    const std::vector<std::string> &files, const bcg::ShearLayout &layout, Config &cfg,
    SelfCalibrationDiagnostics &diagnostics, const bcg::ScratchSpace &scratch, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    const auto budget = cfg.memory_bytes();
    const bool need_bar = cfg.cuts.star_bar_mltp >= 0 || cfg.cuts.expo_bar_mltp > 0;
    const auto source_file = scratch.path("quality_sources.bin");
    const auto bar_file = scratch.path("quality_bars.bin");
    {
        bcg::BinaryWriter<bcg::Source> sources(source_file, 1024 * 1024);
        std::unique_ptr<bcg::BinaryWriter<bcg::BarDatum>> bars;
        if (need_bar) bars.reset(new bcg::BinaryWriter<bcg::BarDatum>(bar_file, 1024 * 1024));
        bcg::for_each_local_source(files, cfg.cuts, layout, comm,
            [&](const bcg::Source &s, const bcg::BarDatum &b) {
                auto raw = s;
                if (cfg.fd_corrected) bcg::restore_raw_shear(raw);
                sources.append(raw);
                if (bars) bars->append(b);
            });
        sources.close(); if (bars) bars->close();
    }
    bcg::StarBarCuts fitted;
    if (need_bar) {
        if (!rank) std::cout << "Fitting stellar bar before pairing" << std::endl;
        fitted = bcg::fit_star_bar(bcg::BarRows(bar_file), cfg.cuts, comm, budget);
    }
    const auto metadata = static_cast<std::uint64_t>(cfg.n_r_bins) *
        (lenses.size() + static_cast<std::uint64_t>(cfg.n_jack) * 8 + 16384) +
        static_cast<std::uint64_t>(lenses.size()) * 128;
    if (metadata > budget / 8)
        throw std::runtime_error("Selected-lens/bin metadata exceeds run.memory-mb budget");
    std::vector<std::unique_ptr<bcg::PairStore>> bins;
    for (int i = 0; i < cfg.n_r_bins; ++i)
        bins.emplace_back(new bcg::PairStore(scratch.path("pairs_" + std::to_string(i) + ".bin"),
            budget / 8, budget / 8 / static_cast<std::size_t>(cfg.n_r_bins), lenses.size(), cfg.n_jack));
    const auto edges = radial_edges(cfg);
    bcg::Cosmology cosmo(cfg.omega_m, cfg.h);
    struct Geometry { double chi = 0, chord = 0; std::array<double, 3> xyz{}; };
    std::vector<Geometry> geometry(lenses.size());
    for (std::size_t i = 0; i < lenses.size(); ++i) {
        auto &g = geometry[i]; const auto &l = lenses[i];
        g.chi = cosmo.chi(l.z);
        g.xyz = bcg::sky_xyz(l.ra, l.dec);
        if (g.chi > 0) g.chord = 2 * std::sin(cfg.rmax / g.chi / 2);
    }
    const std::size_t capacity = std::max<std::size_t>(1, budget / 8 / sizeof(bcg::Source));
    if (capacity > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Source chunk exceeds KD-tree index limit; reduce run.memory-mb");
    std::vector<bcg::Source> chunk; chunk.reserve(capacity);
    std::int64_t local_sources = 0;
    std::unique_ptr<SelfCalibration> selfcal;
    if (cfg.self_calibrate) selfcal.reset(new SelfCalibration(cfg, scratch));
    bool collect_fd = false;
    // ==========================================
    // Function: pair_chunk
    // Method: Use great-circle radii and source-local bearings with the existing response formula.
    // ==========================================
    const auto pair_chunk = [&] {
        if (chunk.empty()) return;
        std::vector<unsigned char> selected(collect_fd ? chunk.size() : 0, 0);
        bcg::KDTree tree(chunk);
        std::vector<int> candidates; candidates.reserve(chunk.size());
        for (std::size_t li = 0; li < lenses.size(); ++li) {
            const auto &lens = lenses[li]; const auto &geo = geometry[li];
            if (geo.chi <= 0) continue;
            tree.radius(geo.xyz, geo.chord, candidates);
            for (int idx : candidates) {
                const auto &s = chunk[idx];
                if (!(s.zp > lens.z + cfg.dzmin)) continue;
                const auto &l = geo.xyz;
                const double xy = s.xyz[0] * l[0] + s.xyz[1] * l[1];
                const double dot = xy + s.xyz[2] * l[2];
                const double cx = s.xyz[1] * l[2] - s.xyz[2] * l[1];
                const double cy = s.xyz[2] * l[0] - s.xyz[0] * l[2];
                const double cz = s.xyz[0] * l[1] - s.xyz[1] * l[0];
                const double radius = geo.chi * std::atan2(std::hypot(cx, cy, cz), dot);
                if (!(radius >= cfg.rmin && radius < cfg.rmax)) continue;
                if (collect_fd) {
                    selected[static_cast<std::size_t>(idx)] = 1;
                    continue;
                }
                const auto it = std::upper_bound(edges.begin(), edges.end(), radius);
                const auto bin = static_cast<std::size_t>(it - edges.begin() - 1);
                // Lens direction in the source's east/north frame; both components
                // carry the same positive cos(source Dec) factor, which cancels.
                const double east = cz;
                const double north = (s.xyz[0] * s.xyz[0] + s.xyz[1] * s.xyz[1]) * l[2]
                    - s.xyz[2] * xy;
                const double phi = std::atan2(north, east);
                const double c2 = std::cos(2 * phi), s2 = std::sin(2 * phi);
                const double c4 = std::cos(4 * phi), s4 = std::sin(4 * phi);
                const double gt = -(s.g1 * c2 + s.g2 * s2);
                const double gx = s.g1 * s2 - s.g2 * c2;
                const double ut = -(s.h1 * c4 + s.h2 * s4);
                const double beta = cosmo.beta(lens.z, s.zp, geo.chi, s.chi);
                bins[bin]->append({gt, gx, beta * (s.de + ut), beta * (s.de - ut),
                                   static_cast<std::int32_t>(li), 0}, radius, labels[li]);
            }
        }
        if (collect_fd)
            for (std::size_t i = 0; i < chunk.size(); ++i)
                if (selected[i]) selfcal->append(chunk[i]);
        chunk.clear();
    };
    // ==========================================
    // Function: scan_sources
    // Method: Read unchanged post-bar rows; apply calibration only during the final ESD pass.
    // ==========================================
    const auto scan_sources = [&] {
        local_sources = 0;
        bcg::BinaryReader<bcg::Source> sources(source_file, 1024 * 1024);
        std::unique_ptr<bcg::BinaryReader<bcg::BarDatum>> bars;
        if (need_bar) bars.reset(new bcg::BinaryReader<bcg::BarDatum>(bar_file, 1024 * 1024));
        bcg::Source source; bcg::BarDatum bar;
        while (sources.next(source)) {
            if (bars && !bars->next(bar)) throw std::runtime_error("Source/bar spool lengths differ");
            if (bars && !fitted.accepts(bar)) continue;
            if (!collect_fd) cfg.calibration.apply(source);
            source.chi = cosmo.chi(source.zp);
            const auto xyz = bcg::sky_xyz(source.ra, source.dec);
            std::copy(xyz.begin(), xyz.end(), source.xyz);
            chunk.push_back(source); ++local_sources;
            if (chunk.size() == capacity) pair_chunk();
        }
        if (bars && bars->next(bar)) throw std::runtime_error("Source/bar spool lengths differ");
        pair_chunk();
    };
    if (selfcal) {
        if (!rank) std::cout << "Collecting raw all-R source union for Cartesian FD self-calibration" << std::endl;
        collect_fd = true;
        scan_sources();
        cfg.calibration = selfcal->finish(comm, diagnostics);
        cfg.self_calibration_sources = selfcal->selected_sources();
        selfcal.reset();
        collect_fd = false;
    }
    if (!rank) std::cout << "Writing bounded pair spools" << std::endl;
    scan_sources();
    std::vector<bcg::Source>().swap(chunk);
    bcg::fs::remove(source_file); if (need_bar) bcg::fs::remove(bar_file);
    std::int64_t local_pairs = 0;
    for (auto &bin : bins) { bin->finish(); local_pairs += bin->count; }
    std::int64_t total_sources = 0, total_pairs = 0;
    MPI_Allreduce(&local_sources, &total_sources, 1, MPI_INT64_T, MPI_SUM, comm);
    MPI_Allreduce(&local_pairs, &total_pairs, 1, MPI_INT64_T, MPI_SUM, comm);
    if (!rank) std::cout << "Sources after cuts: " << total_sources << '\n'
                         << "BCG-source pairs: " << total_pairs << std::endl;
    if (!total_sources) throw std::runtime_error("No sources survive quality cuts");
    if (!total_pairs) throw std::runtime_error("No lens-source pairs");
    return bins;
}
} // namespace bcg_esd
