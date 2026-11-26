#include "bcg_esd.hpp"
#include "bcg_star_bar.hpp"
#include "bcg_physics.hpp"
#include "selfcalibration.hpp"
#include <iostream>

namespace bcg_esd {
namespace {
using bcg::LensGeometry;

// A source can be in both search apertures only if their chord-distance
// bounds intersect. Give up the test early in crowded fields.
bool apertures_are_disjoint(const std::vector<LensGeometry> &geometry) {
    std::vector<const LensGeometry *> sorted;
    sorted.reserve(geometry.size());
    double max_chord = 0;
    for (const auto &lens : geometry) {
        if (lens.chi <= 0) continue;
        if (!std::isfinite(lens.chord)) return false;
        sorted.push_back(&lens);
        max_chord = std::max(max_chord, std::abs(lens.chord));
    }
    std::sort(sorted.begin(), sorted.end(), [](const auto *a, const auto *b) { return a->xyz[0] < b->xyz[0]; });
    std::size_t remaining_checks = 8 * sorted.size();
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        const auto &a = *sorted[i];
        const double radius = std::abs(a.chord) + 1e-12;
        for (std::size_t j = i + 1; j < sorted.size(); ++j) {
            const auto &b = *sorted[j];
            const double dx = b.xyz[0] - a.xyz[0];
            if (dx > radius + max_chord) break;
            if (!remaining_checks) return false;
            --remaining_checks;
            const double dy = b.xyz[1] - a.xyz[1], dz = b.xyz[2] - a.xyz[2];
            const double sum = radius + std::abs(b.chord);
            if (dx * dx + dy * dy + dz * dz <= sum * sum) return false;
        }
    }
    return true;
}
} // namespace

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
    const std::vector<bcg::Lens> &lenses,
    const std::vector<std::string> &files, const bcg::ShearLayout &layout, Config &cfg,
    const bcg::ScratchSpace &scratch, std::int64_t &sources_after_cuts, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    const auto budget = cfg.memory_bytes();
    const bool need_bar = cfg.cuts.bar_enabled();
    const std::size_t capacity = std::max<std::size_t>(1, budget / 8 / sizeof(bcg::Source));
    if (capacity > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Source chunk exceeds KD-tree index limit; reduce run.memory-mb");
    const auto source_file = [&](std::size_t part) {
        return scratch.path("quality_sources_" + std::to_string(part) + ".bin");
    };
    const auto bar_file = scratch.path("quality_bars.bin");
    const auto restore = [&](bcg::Source source) {
        if (cfg.fd_corrected) bcg::restore_raw_shear(source);
        return source;
    };
    constexpr std::size_t io_bytes = 1024 * 1024;
    bcg::StarBarCuts fitted;
    std::size_t quality_parts = 0;
    if (need_bar) {
        std::unique_ptr<bcg::BinaryWriter<bcg::SourceRecord>> sources;
        bcg::BinaryWriter<bcg::BarDatum> bars(bar_file, io_bytes);
        std::size_t remaining = 0;
        bcg::for_each_local_source(files, cfg.cuts, layout, comm,
            [&](const bcg::Source &s, const bcg::BarDatum &b) {
                // Bound the overlap of the input cache and new pair files.
                if (!remaining) {
                    if (sources) sources->close();
                    sources = std::make_unique<bcg::BinaryWriter<bcg::SourceRecord>>(source_file(quality_parts++), io_bytes);
                    remaining = std::max<std::size_t>(1, capacity / 4);
                }
                sources->append(restore(s));
                bars.append(b);
                --remaining;
            });
        if (sources) sources->close();
        bars.close();
        if (!rank) std::cout << "Fitting stellar bar before pairing" << std::endl;
        fitted = bcg::fit_star_bar(bcg::BarRows(bar_file), cfg.cuts, comm, budget);
    }
    const auto metadata = static_cast<std::uint64_t>(cfg.n_r_bins) *
        (bcg::lens_mask_words(lenses.size()) * sizeof(std::uint64_t) + 16384) +
        static_cast<std::uint64_t>(lenses.size()) * 128;
    if (metadata > budget / 8)
        throw std::runtime_error("Selected-lens/bin metadata exceeds run.memory-mb budget");
    const auto edges = radial_edges(cfg);
    bcg::Cosmology cosmo(cfg.cosmology.omega_m, cfg.cosmology.h);
    std::vector<LensGeometry> geometry(lenses.size());
    for (std::size_t i = 0; i < lenses.size(); ++i) {
        auto &g = geometry[i]; const auto &l = lenses[i];
        g.z = l.z;
        g.chi = cosmo.chi(l.z);
        g.xyz = bcg::sky_xyz(l.ra, l.dec);
        if (g.chi > 0) g.chord = 2 * std::sin(cfg.rmax / g.chi / 2);
    }
    std::shared_ptr<bcg::SourceMomentStore> matched;
    if (cfg.self_calibrate || !apertures_are_disjoint(geometry))
        matched = std::make_shared<bcg::SourceMomentStore>(scratch.path("matched_sources.bin"),
                                                           io_bytes, cfg.self_calibrate, cosmo);
    const auto bin_write_bytes = std::min(io_bytes / 4, budget / 8 / static_cast<std::size_t>(cfg.n_r_bins));
    std::vector<std::unique_ptr<bcg::PairStore>> bins;
    for (int i = 0; i < cfg.n_r_bins; ++i)
        bins.emplace_back(new bcg::PairStore(scratch.path("pairs_" + std::to_string(i) + ".bin"),
            io_bytes, bin_write_bytes, lenses.size(), 1, false, matched));
    std::vector<bcg::Source> chunk; chunk.reserve(capacity);
    std::int64_t local_sources = 0;
    std::unique_ptr<SelfCalibration> selfcal;
    if (cfg.self_calibrate) selfcal.reset(new SelfCalibration(cfg, scratch));
    std::vector<std::uint32_t> preorder;
    std::vector<unsigned char> selected;
    // ==========================================
    // Function: pair_chunk
    // Method: Use great-circle radii and source-local bearings with the existing response formula.
    // ==========================================
    const auto pair_chunk = [&] {
        if (chunk.empty()) return;
        const auto block = matched ? matched->block_id() : 0;
        bcg::KDTree tree(chunk);
        if (matched) {
            selected.assign(chunk.size(), 0);
            preorder.resize(chunk.size());
            for (std::size_t i = 0; i < chunk.size(); ++i)
                preorder[tree.preorder_point(i)] = static_cast<std::uint32_t>(i);
        }
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
                const auto it = std::upper_bound(edges.begin(), edges.end(), radius);
                const auto bin = static_cast<std::size_t>(it - edges.begin() - 1);
                const auto lens_index = static_cast<std::int32_t>(li);
                if (matched) {
                    selected[idx] = 1;
                    bins[bin]->append_link(block, preorder[idx], lens_index, radius);
                } else {
                    const double phi = bcg::source_bearing(s.xyz, l, xy, cz);
                    const double beta = bcg::Cosmology::beta(lens.z, s.zp, geo.chi, s.chi);
                    bins[bin]->append(bcg::project_pair({s.g1, s.g2, s.de, s.h1, s.h2},
                                                      phi, beta, lens_index), radius);
                }
            }
        }
        // Preserve the original source order for the FD-bin sums.
        if (selfcal)
            for (std::size_t i = 0; i < chunk.size(); ++i)
                if (selected[i]) selfcal->append(chunk[i]);
        if (matched) {
            for (std::size_t i = 0; i < chunk.size(); ++i) {
                const auto idx = tree.preorder_point(i);
                if (selected[idx])
                    matched->append(chunk[idx], static_cast<std::uint32_t>(i));
            }
            matched->finish_block(chunk.size());
        }
        chunk.clear();
    };
    const auto accept = [&](bcg::Source source) {
        if (!selfcal) cfg.calibration.apply(source);
        source.chi = cosmo.chi(source.zp);
        const auto xyz = bcg::sky_xyz(source.ra, source.dec);
        std::copy(xyz.begin(), xyz.end(), source.xyz);
        chunk.push_back(source); ++local_sources;
        if (chunk.size() == capacity) pair_chunk();
    };
    if (!rank) std::cout << (selfcal ? "Matching sources and collecting the raw all-R FD union"
                                    : "Writing bounded pair spools") << std::endl;
    if (need_bar) {
        bcg::BinaryReader<bcg::BarDatum> bars(bar_file, io_bytes);
        bcg::BarDatum bar;
        for (std::size_t part = 0; part < quality_parts; ++part) {
            const auto path = source_file(part);
            {
                bcg::BinaryReader<bcg::SourceRecord> sources(path, io_bytes);
                bcg::SourceRecord record;
                while (sources.next(record)) {
                    if (!bars.next(bar)) throw std::runtime_error("Source/bar spool lengths differ");
                    if (fitted.accepts(bar)) accept(bcg::Source(record));
                }
            }
            bcg::fs::remove(path);
        }
        if (bars.next(bar)) throw std::runtime_error("Source/bar spool lengths differ");
    } else {
        bcg::for_each_local_source(files, cfg.cuts, layout, comm,
            [&](const bcg::Source &source, const bcg::BarDatum &) { accept(restore(source)); });
    }
    if (need_bar) bcg::fs::remove(bar_file);
    pair_chunk();
    std::vector<bcg::Source>().swap(chunk);
    std::vector<std::uint32_t>().swap(preorder);
    std::vector<unsigned char>().swap(selected);
    std::int64_t local_pairs = 0;
    for (auto &bin : bins) { bin->finish(); local_pairs += bin->count; }
    if (matched) matched->finish(std::move(geometry));
    if (selfcal) {
        cfg.calibration = selfcal->finish(comm);
        cfg.self_calibration_sources = selfcal->selected_sources();
        matched->calibrate(cfg.calibration);
    }
    std::int64_t total_pairs = 0;
    MPI_Allreduce(&local_sources, &sources_after_cuts, 1, MPI_INT64_T, MPI_SUM, comm);
    MPI_Allreduce(&local_pairs, &total_pairs, 1, MPI_INT64_T, MPI_SUM, comm);
    if (!rank) std::cout << "Sources after cuts: " << sources_after_cuts << '\n'
                         << "Lens-source pairs: " << total_pairs << std::endl;
    if (!sources_after_cuts) throw std::runtime_error("No sources survive quality cuts");
    if (!total_pairs) throw std::runtime_error("No lens-source pairs");
    return bins;
}
} // namespace bcg_esd
