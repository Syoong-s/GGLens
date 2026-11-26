#include "bcg_esd.hpp"
#include "bcg_mpi_rows.hpp"
#include "bcg_physics.hpp"
#include <iostream>
#include <random>

namespace bcg_esd {
namespace {
using Point = std::array<double, 3>;

// ==========================================
// Function: dot
// Method: Compare spherical distances without computing an angle.
// ==========================================
double dot(const Point &a, const Point &b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// ==========================================
// Function: seed_centers
// Method: Initialize spherical k-means on a deterministic sample of the merged catalog.
// ==========================================
std::vector<Point> seed_centers(const std::vector<Point> &points, int regions) {
    std::mt19937 rng(1907);
    std::vector<Point> sample;
    const auto sample_size = std::min(points.size(), std::max<std::size_t>(8192, regions));
    std::sample(points.begin(), points.end(), std::back_inserter(sample), sample_size, rng);
    std::vector<double> distance(sample.size(), 2);
    std::vector<Point> centers;
    centers.push_back(sample[std::uniform_int_distribution<std::size_t>(0, sample.size() - 1)(rng)]);
    while (static_cast<int>(centers.size()) < regions) {
        double sum = 0;
        for (std::size_t i = 0; i < sample.size(); ++i) {
            distance[i] = std::min(distance[i], std::max(0.0, 1.0 - dot(sample[i], centers.back())));
            sum += distance[i];
        }
        if (!(sum > 0)) break;
        double step = std::uniform_real_distribution<double>(0.0, sum)(rng);
        auto pos = sample.size() - 1;
        for (std::size_t i = 0; i < sample.size(); ++i) {
            step -= distance[i];
            if (step < 0) { pos = i; break; }
        }
        centers.push_back(sample[pos]);
    }
    return centers;
}

// ==========================================
// Function: spherical_regions
// Method: Assign rank-local sky positions, reduce centroid sums, and compact empty cells.
// ==========================================
std::vector<int> spherical_regions(const std::vector<Point> &points, int regions, MPI_Comm comm) {
    int rank = 0, ranks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    auto centers = rank == 0 ? seed_centers(points, regions) : std::vector<Point>{};
    regions = static_cast<int>(centers.size());
    MPI_Bcast(&regions, 1, MPI_INT, 0, comm);
    centers.resize(regions);
    MPI_Bcast(centers.data(), regions * sizeof(Point), MPI_BYTE, 0, comm);

    const auto [begin, end] = bcg::partition(points.size(), rank, ranks);
    std::vector<int> local(end - begin, -1);
    std::vector<double> sums(4 * regions);
    for (int iteration = 0; iteration < 32; ++iteration) {
        int changed = 0;
        std::fill(sums.begin(), sums.end(), 0);
        for (auto i = begin; i < end; ++i) {
            const Point &point = points[i];
            int best = 0;
            double similarity = dot(point, centers.front());
            for (int r = 1; r < regions; ++r) {
                const double current = dot(point, centers[r]);
                if (current > similarity) { similarity = current; best = r; }
            }
            changed |= local[i - begin] != best;
            local[i - begin] = best;
            for (int j = 0; j < 3; ++j) sums[4*best + j] += point[j];
            ++sums[4*best + 3];
        }
        MPI_Allreduce(MPI_IN_PLACE, sums.data(), 4 * regions, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(MPI_IN_PLACE, &changed, 1, MPI_INT, MPI_MAX, comm);
        if (!changed) break;
        for (int r = 0; r < regions; ++r) {
            const auto *sum = sums.data() + 4*r;
            const double norm = std::hypot(sum[0], sum[1], sum[2]);
            if (norm > 0) centers[r] = {sum[0]/norm, sum[1]/norm, sum[2]/norm};
        }
    }
    std::vector<int> compact(regions, -1);
    int next = 0;
    for (int r = 0; r < regions; ++r)
        if (sums[4*r + 3] > 0) compact[r] = next++;
    for (auto &label : local) label = compact[label];
    return bcg::allgather_rows(local, comm);
}
} // namespace

// ==========================================
// Function: make_lens_stack
// Method: Use the same participating group sample to set both sky regions and the fixed model inputs.
// ==========================================
LensStack make_lens_stack(const std::vector<bcg::Lens> &lenses,
                          const std::vector<std::unique_ptr<bcg::PairStore>> &bins,
                          const Config &cfg, MPI_Comm comm) {
    std::vector<std::uint64_t> active(bcg::lens_mask_words(lenses.size()), 0);
    for (const auto &bin : bins) {
        auto count = bin->count;
        MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT64_T, MPI_SUM, comm);
        if (count >= cfg.pdf_bins + 2LL)
            for (std::size_t i = 0; i < active.size(); ++i) active[i] |= bin->present[i];
        std::vector<std::uint64_t>().swap(bin->present);
    }
    MPI_Allreduce(MPI_IN_PLACE, active.data(), static_cast<int>(active.size()),
                  MPI_UINT64_T, MPI_BOR, comm);
    const auto has_pairs = [&](std::size_t i) {
        return (active[i / bcg::lenses_per_word] >> (i % bcg::lenses_per_word)) & std::uint64_t{1};
    };
    LensStack stack;
    auto &sample = stack.sample;
    stack.labels.assign(lenses.size(), -1);
    std::vector<Point> points;
    std::size_t reliable = 0;
    const bool report_photometry = cfg.center == bcg::Center::bcg && cfg.lens_columns.galaxy.photometry_select != 0;
    for (std::size_t i = 0; i < lenses.size(); ++i) {
        if (!has_pairs(i)) continue;
        const auto &lens = lenses[i];
        points.push_back(bcg::sky_xyz(lens.ra, lens.dec));
        sample.redshift += lens.z;
        if (cfg.center == bcg::Center::bcg) {
            sample.stellar_mass += lens.stellar_mass;
            if (report_photometry) reliable += lens.photometry_select >= 1;
        }
    }
    sample.n_lenses = points.size();
    if (sample.n_lenses < 2) throw std::runtime_error("Jackknife needs at least two groups with usable pairs");
    sample.stellar_mass *= cfg.cosmology.h / sample.n_lenses; // Convert physical M_sun to h^-1 M_sun.
    sample.redshift /= sample.n_lenses;

    const int requested = std::min(cfg.n_jack, static_cast<int>(sample.n_lenses));
    const auto labels = spherical_regions(points, requested, comm);
    sample.regions = 1 + *std::max_element(labels.begin(), labels.end());
    if (sample.regions < 2) throw std::runtime_error("Lens sky positions cannot form jackknife regions");
    std::vector<int> counts(sample.regions);
    std::size_t pos = 0;
    for (std::size_t i = 0; i < lenses.size(); ++i) {
        if (!has_pairs(i)) continue;
        stack.labels[i] = labels[pos++];
        ++counts[stack.labels[i]];
    }
    const auto bounds = std::minmax_element(counts.begin(), counts.end());
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (!rank) {
        std::cout << "Combined " << bcg::center_name(cfg.center) << " stack: " << sample.n_lenses
                  << " groups, mean z=" << sample.redshift;
        if (report_photometry) std::cout << ", select>=1: " << reliable;
        if (cfg.center == bcg::Center::bcg)
            std::cout << ", mean stellar mass [M_sun/h]=" << sample.stellar_mass;
        std::cout << "\nJackknife regions: " << sample.regions << '/' << cfg.n_jack
                  << ", groups per region " << *bounds.first << " to " << *bounds.second << std::endl;
    }
    stack.photometry_selected = report_photometry ? static_cast<std::int64_t>(reliable) : -1;
    stack.groups_per_region = std::move(counts);
    return stack;
}
} // namespace bcg_esd
