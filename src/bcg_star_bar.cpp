#include "bcg_star_bar.hpp"

#include <mpi.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bcg {
namespace detail {
constexpr int size_bins = 100;
constexpr int mag_bins = 20;
constexpr double size_min = -2.0;
constexpr double size_width = 0.04;
constexpr double mag_min = 10.0;
constexpr double mag_width = 1.0;

// ==========================================
// Function: mag_bin
// Method: Convert the reference's truncating magnitude bins to zero-based indices.
// ==========================================
int mag_bin(double mag) {
    return static_cast<int>((mag - mag_min) / mag_width);
}

// ==========================================
// Function: size_bin
// Method: Return the half-open size histogram bin, or -1 when out of range.
// ==========================================
int size_bin(double size) {
    if (size < size_min || size >= size_min + size_width * size_bins) return -1;
    return static_cast<int>((size - size_min) / size_width);
}

// ==========================================
// Function: checked_mpi_count
// Method: Keep dynamic reductions within the MPI int-count API.
// ==========================================
int checked_mpi_count(std::size_t n) {
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Star-bar MPI reduction exceeds int count");
    return static_cast<int>(n);
}

// ==========================================
// Function: reduce_counts
// Method: Sum a same-sized integer vector on all ranks.
// ==========================================
void reduce_counts(std::vector<std::int64_t> &values, MPI_Comm comm) {
    MPI_Allreduce(MPI_IN_PLACE, values.data(), checked_mpi_count(values.size()),
                  MPI_INT64_T, MPI_SUM, comm);
}

// ==========================================
// Function: reduce_doubles
// Method: Sum a same-sized double vector on all ranks.
// ==========================================
void reduce_doubles(std::vector<double> &values, MPI_Comm comm) {
    MPI_Allreduce(MPI_IN_PLACE, values.data(), checked_mpi_count(values.size()),
                  MPI_DOUBLE, MPI_SUM, comm);
}

// ==========================================
// Function: size_center
// Method: Convert a zero-based histogram index to its bin center.
// ==========================================
double size_center(int i) {
    return size_min + (static_cast<double>(i) + 0.5) * size_width;
}

// ==========================================
// Function: peak_concentration
// Method: Count the peak size bin and its two neighbors per magnitude bin.
// ==========================================
double peak_concentration(const std::int64_t *hist, std::int64_t count, int peak) {
    std::int64_t three = hist[peak];
    if (peak > 0) three += hist[peak - 1];
    if (peak + 1 < size_bins) three += hist[peak + 1];
    return count > 0 ? static_cast<double>(three) / count : 0.0;
}
} // namespace detail

// ==========================================
// Function: fit_global_star_cut
// Method: Reproduce the reference's global size-magnitude star peak and sigma cut.
// ==========================================
double fit_global_star_cut(const BarRows &bar,
                                  double multiplier, MPI_Comm comm) {
    using namespace detail;
    std::vector<std::int64_t> hist(size_bins * mag_bins, 0), mag_count(mag_bins, 0);
    bar.each([&](const BarDatum &b) {
        const int j = mag_bin(b.mag_i);
        if (j < 0 || j >= mag_bins) return;
        ++mag_count[j];
        const int i = size_bin(b.size_rel);
        if (i >= 0) ++hist[j * size_bins + i];
    });
    reduce_counts(hist, comm);
    reduce_counts(mag_count, comm);

    constexpr int limit = static_cast<int>((-0.5 - size_min) / size_width) + 1;
    double max_concentration = 0, initial = 0.5;
    int best_mag = -1;
    for (int j = 0; j < mag_bins; ++j) {
        if (mag_count[j] < 100) continue;
        const std::int64_t *h = hist.data() + j * size_bins;
        int peak = limit;
        for (int i = limit + 1; i < size_bins; ++i)
            if (h[i] > h[peak]) peak = i;
        const double concentration = peak_concentration(h, mag_count[j], peak);
        if (concentration > max_concentration) {
            max_concentration = concentration;
            best_mag = j;
            initial = size_center(peak);
        }
    }
    std::array<unsigned char, mag_bins> active{};
    if (best_mag >= 0)
        for (int j = 0; j < mag_bins; ++j) {
            if (mag_count[j] < 100) continue;
            const std::int64_t *h = hist.data() + j * size_bins;
            int peak = limit;
            for (int i = limit + 1; i < size_bins; ++i)
                if (h[i] > h[peak]) peak = i;
            if (std::abs(size_center(peak) - initial) <= 0.05 &&
                peak_concentration(h, mag_count[j], peak) >= 0.6)
                active[j] = 1;
        }

    double local_sum = 0;
    std::int64_t local_count = 0;
    bar.each([&](const BarDatum &b) {
        const int j = mag_bin(b.mag_i);
        if (j >= 0 && j < mag_bins && active[j] &&
            std::abs(b.size_rel - initial) < 0.1) {
            local_sum += b.size_rel;
            ++local_count;
        }
    });
    double global_sum = 0;
    std::int64_t global_count = 0;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT64_T, MPI_SUM, comm);
    if (global_count <= 1)
        return best_mag >= 0 ? initial + 0.05 * multiplier : 0.6;

    const double mean = global_sum / static_cast<double>(global_count);
    double local_sq = 0;
    bar.each([&](const BarDatum &b) {
        const int j = mag_bin(b.mag_i);
        if (j >= 0 && j < mag_bins && active[j] &&
            std::abs(b.size_rel - initial) < 0.1)
            local_sq += (b.size_rel - mean) * (b.size_rel - mean);
    });
    double global_sq = 0;
    MPI_Allreduce(&local_sq, &global_sq, 1, MPI_DOUBLE, MPI_SUM, comm);
    return mean + multiplier * std::sqrt(global_sq / (global_count - 1));
}

// ==========================================
// Function: global_exposure_ids
// Method: Collect unique physical exposure IDs from all MPI ranks.
// ==========================================
std::vector<int> global_exposure_ids(const BarRows &bar, MPI_Comm comm) {
    constexpr std::size_t words = (10000001 + 63) / 64;
    std::vector<std::uint64_t> bits(words, 0);
    bar.each([&](const BarDatum &b) {
        if (b.exposure_id < 1 || b.exposure_id > 10000000)
            throw std::runtime_error("Invalid exposure ID in star-bar metadata");
        bits[static_cast<std::size_t>(b.exposure_id) / 64] |= std::uint64_t{1} << (b.exposure_id % 64);
    });
    MPI_Allreduce(MPI_IN_PLACE, bits.data(), static_cast<int>(words), MPI_UINT64_T, MPI_BOR, comm);
    std::vector<int> ids;
    std::size_t n = 0;
    for (auto word : bits) n += static_cast<std::size_t>(__builtin_popcountll(word));
    ids.reserve(n);
    for (int id = 1; id <= 10000000; ++id)
        if (bits[static_cast<std::size_t>(id) / 64] & (std::uint64_t{1} << (id % 64))) ids.push_back(id);
    return ids;
}

// ==========================================
// Structure: ExposureFits
// Method: Keep batch-local means and scatter before applying the global fallback.
// ==========================================
struct ExposureFits {
    std::vector<double> means, sigmas;
    std::vector<unsigned char> skipped;
};

// ==========================================
// Function: fit_exposure_star_cuts
// Method: Fit the reference's per-exposure stellar peak, clipping, and fallback.
// ==========================================
ExposureFits fit_exposure_star_cuts(
    const BarRows &bar, const std::vector<int> &ids,
    MPI_Comm comm) {
    using namespace detail;
    const std::size_t n = ids.size();
    if (!n) return {};
    checked_mpi_count(n * size_bins * mag_bins);
    std::unordered_map<int, std::size_t> slot;
    slot.reserve(n);
    for (std::size_t e = 0; e < n; ++e) slot.emplace(ids[e], e);
    std::vector<std::int64_t> hist(n * size_bins * mag_bins, 0), mag_count(n * mag_bins, 0);
    bar.each([&](const BarDatum &b) {
        if (b.snr <= 40.0) return;
        const std::size_t e = slot.at(b.exposure_id);
        const int j = mag_bin(b.mag_i), i = size_bin(b.size_rel);
        if (j < 0 || j >= mag_bins) return;
        ++mag_count[e * mag_bins + j];
        if (i >= 0) ++hist[(e * mag_bins + j) * size_bins + i];
    });
    reduce_counts(hist, comm);
    reduce_counts(mag_count, comm);

    constexpr int star_first = static_cast<int>((-0.5 - size_min) / size_width);
    constexpr int star_last = static_cast<int>((0.2 - size_min) / size_width);
    std::vector<double> initial(n, 0.5), mean(n, 0.5), stddev(n, 0.05);
    std::vector<unsigned char> active(n * mag_bins, 0), fallback(n, 0), skip(n, 0);
    for (std::size_t e = 0; e < n; ++e) {
        double best_concentration = 0;
        int best_mag = -1;
        for (int j = 0; j < mag_bins; ++j) {
            const std::int64_t count = mag_count[e * mag_bins + j];
            if (count < 100) continue;
            const std::int64_t *h = hist.data() + (e * mag_bins + j) * size_bins;
            int peak = star_first;
            for (int i = star_first + 1; i <= star_last; ++i)
                if (h[i] > h[peak]) peak = i;
            const double concentration = peak_concentration(h, count, peak);
            if (concentration > best_concentration) {
                best_concentration = concentration;
                best_mag = j;
                initial[e] = size_center(peak);
            }
        }
        if (best_mag < 0) continue;
        for (int j = 0; j < mag_bins; ++j) {
            const std::int64_t count = mag_count[e * mag_bins + j];
            if (count < 100) continue;
            const std::int64_t *h = hist.data() + (e * mag_bins + j) * size_bins;
            int peak = star_first;
            for (int i = star_first + 1; i <= star_last; ++i)
                if (h[i] > h[peak]) peak = i;
            if (std::abs(size_center(peak) - initial[e]) <= 0.05 &&
                peak_concentration(h, count, peak) >= 0.6)
                active[e * mag_bins + j] = 1;
        }
    }
    hist.clear();
    hist.shrink_to_fit();
    mag_count.clear();
    mag_count.shrink_to_fit();

    std::vector<double> sum(n, 0), square(n, 0), limit(n, 0);
    std::vector<std::int64_t> count(n, 0);
    bar.each([&](const BarDatum &b) {
        if (b.snr <= 40.0) return;
        const std::size_t e = slot.at(b.exposure_id);
        const int j = mag_bin(b.mag_i);
        if (j < 0 || j >= mag_bins || !active[e * mag_bins + j] ||
            std::abs(b.size_rel - initial[e]) >= 0.1) return;
        sum[e] += b.size_rel;
        ++count[e];
    });
    reduce_doubles(sum, comm);
    reduce_counts(count, comm);
    for (std::size_t e = 0; e < n; ++e) {
        fallback[e] = count[e] <= 1;
        skip[e] = fallback[e];
        if (!fallback[e]) {
            mean[e] = sum[e] / count[e];
            stddev[e] = 0.1;
        }
    }

    std::fill(sum.begin(), sum.end(), 0);
    std::fill(count.begin(), count.end(), 0);
    bar.each([&](const BarDatum &b) {
        if (b.snr <= 40.0) return;
        const std::size_t e = slot.at(b.exposure_id);
        const int j = mag_bin(b.mag_i);
        if (!fallback[e] || j < 0 || j >= mag_bins ||
            std::abs(b.size_rel - initial[e]) >= 0.15) return;
        sum[e] += b.size_rel;
        ++count[e];
    });
    reduce_doubles(sum, comm);
    reduce_counts(count, comm);
    for (std::size_t e = 0; e < n; ++e)
        if (fallback[e] && count[e] > 1) {
            mean[e] = sum[e] / count[e];
            stddev[e] = 0.15;
            skip[e] = 0;
        }

    for (int iteration = 1; iteration <= 3; ++iteration) {
        for (std::size_t e = 0; e < n; ++e) {
            if (skip[e]) continue;
            const double window = fallback[e] ? 0.15 : 0.1;
            if (iteration == 1) limit[e] = window;
            else if (iteration == 2) limit[e] = 0.5 * window;
            else limit[e] = std::max(std::min(3.0 * stddev[e], 0.5 * window), 0.015);
        }
        std::fill(sum.begin(), sum.end(), 0);
        std::fill(count.begin(), count.end(), 0);
        bar.each([&](const BarDatum &b) {
            if (b.snr <= 40.0) return;
            const std::size_t e = slot.at(b.exposure_id);
            if (skip[e]) return;
            const int j = mag_bin(b.mag_i);
            if (j < 0 || j >= mag_bins ||
                (!fallback[e] && !active[e * mag_bins + j]) ||
                std::abs(b.size_rel - mean[e]) >= limit[e]) return;
            sum[e] += b.size_rel;
            ++count[e];
        });
        reduce_doubles(sum, comm);
        reduce_counts(count, comm);
        for (std::size_t e = 0; e < n; ++e) {
            if (skip[e]) continue;
            if (count[e] > 1) mean[e] = sum[e] / count[e];
            else skip[e] = 1;
        }
        std::fill(square.begin(), square.end(), 0);
        std::fill(count.begin(), count.end(), 0);
        bar.each([&](const BarDatum &b) {
            if (b.snr <= 0.0) return;
            const std::size_t e = slot.at(b.exposure_id);
            if (skip[e]) return;
            const int j = mag_bin(b.mag_i);
            if (j < 0 || j >= mag_bins ||
                (!fallback[e] && !active[e * mag_bins + j]) ||
                std::abs(b.size_rel - mean[e]) >= limit[e]) return;
            const double difference = b.size_rel - mean[e];
            square[e] += difference * difference;
            ++count[e];
        });
        reduce_doubles(square, comm);
        reduce_counts(count, comm);
        for (std::size_t e = 0; e < n; ++e) {
            if (skip[e]) continue;
            stddev[e] = count[e] > 1 ? std::sqrt(square[e] / (count[e] - 1)) : 0.05;
            stddev[e] = std::max(stddev[e], 0.005);
        }
    }
    return {std::move(mean), std::move(stddev), std::move(skip)};
}

// ==========================================
// Function: fit_star_bar
// Method: Fit once globally or batch exposure histograms, retaining the all-exposure fallback.
// ==========================================
StarBarCuts fit_star_bar(const BarRows &bar, const Cuts &cuts, MPI_Comm comm, std::size_t budget) {
    StarBarCuts result;
    if (cuts.expo_bar_mltp <= 0 && cuts.star_bar_mltp < 0) return result;
    result.enabled = true;
    result.exposure = cuts.expo_bar_mltp > 0;
    if (!result.exposure) {
        result.global = fit_global_star_cut(bar, cuts.star_bar_mltp, comm); return result;
    }
    result.ids = global_exposure_ids(bar, comm);
    const std::size_t n = result.ids.size();
    if (n > budget / 8 / (sizeof(int) + sizeof(double) + sizeof(unsigned char)))
        throw std::runtime_error("Exposure-cut metadata exceeds run.memory-mb budget");
    result.thresholds.resize(n);
    std::vector<unsigned char> skipped(n, 0);
    const std::size_t batch_size = std::max<std::size_t>(1, budget / 8 / 17000);
    double mean_sum = 0, sigma_sum = 0;
    std::size_t success = 0;
    for (std::size_t first = 0; first < n; first += batch_size) {
        const auto last = std::min(first + batch_size, n);
        std::vector<int> ids(result.ids.begin() + first, result.ids.begin() + last);
        const auto fits = fit_exposure_star_cuts(bar.subset(ids.front(), ids.back()), ids, comm);
        for (std::size_t i = first; i < last; ++i) {
            const auto j = i - first;
            skipped[i] = fits.skipped[j];
            if (!skipped[i]) {
                mean_sum += fits.means[j]; sigma_sum += fits.sigmas[j]; ++success;
                result.thresholds[i] = fits.means[j] + cuts.expo_bar_mltp * fits.sigmas[j];
            }
        }
    }
    const double fallback_mean = success ? mean_sum / success : 0.5;
    const double fallback_sigma = success ? sigma_sum / success : 0.05;
    for (std::size_t i = 0; i < n; ++i) if (skipped[i])
        result.thresholds[i] = fallback_mean + cuts.expo_bar_mltp * fallback_sigma;
    return result;
}

// ==========================================
// Function: StarBarCuts::accepts
// Method: Apply the original precedence, size inequality and global nonpositive-SNR rule.
// ==========================================
bool StarBarCuts::accepts(const BarDatum &b) const {
    if (!enabled) return true;
    if (!exposure) return !(b.size_rel <= global && b.snr > 0);
    const auto it = std::lower_bound(ids.begin(), ids.end(), b.exposure_id);
    if (it == ids.end() || *it != b.exposure_id)
        throw std::runtime_error("Missing exposure threshold in star-bar fit");
    return !(b.size_rel <= thresholds[static_cast<std::size_t>(it - ids.begin())]);
}

// ==========================================
// Function: apply_star_bar
// Method: Apply the reference's exposure bar when enabled, otherwise global bar.
// ==========================================
void apply_star_bar(LocalSources &local, const Cuts &cuts, MPI_Comm comm) {
    if (cuts.expo_bar_mltp <= 0 && cuts.star_bar_mltp < 0) return;
    if (local.bar_data.size() != local.sources.size())
        throw std::runtime_error("Star-bar metadata does not match source count");
    const auto fitted = fit_star_bar(BarRows(local.bar_data), cuts, comm, 512ULL * 1024 * 1024);
    std::size_t retained = 0;
    for (std::size_t i = 0; i < local.sources.size(); ++i)
        if (fitted.accepts(local.bar_data[i])) local.sources[retained++] = local.sources[i];
    local.sources.resize(retained);
    local.bar_data.clear();
    local.bar_data.shrink_to_fit();
}
} // namespace bcg
