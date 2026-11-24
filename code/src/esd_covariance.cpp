#include "bcg_esd.hpp"
#include <numeric>

namespace bcg_esd {
namespace {
// ==========================================
// Function: replica_mean
// Method: Average the common set of delete-one estimates, propagating unavailable replicas.
// ==========================================
double replica_mean(const std::vector<double> &values) {
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

// ==========================================
// Function: jackknife_entry
// Method: Apply the common-region (K-1)/K normalization to a radial covariance entry.
// ==========================================
double jackknife_entry(const std::vector<double> &a, double a_mean,
                       const std::vector<double> &b, double b_mean) {
    double sum = 0;
    for (std::size_t k = 0; k < a.size(); ++k) sum += (a[k] - a_mean) * (b[k] - b_mean);
    return sum * (a.size() - 1.0) / a.size();
}

} // namespace

// ==========================================
// Function: jackknife_covariance
// Method: Retain every regional deletion for all fitted radii and compute the full matrix.
// ==========================================
Covariance jackknife_covariance(std::vector<Result> &results, int regions, bool diagonal_only) {
    const auto n = results.size();
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    Covariance cov;
    cov.matrix.assign(n*n, nan);
    std::vector<double> means(n, nan);
    for (std::size_t i = 0; i < n; ++i) {
        auto &r = results[i];
        if (r.tangential_replicas.size() != static_cast<std::size_t>(regions) ||
            r.cross_replicas.size() != static_cast<std::size_t>(regions))
            throw std::runtime_error("Inconsistent jackknife region count");
        means[i] = replica_mean(r.tangential_replicas);
        if (std::isfinite(r.measured) && std::isfinite(means[i]))
            cov.fit_bins.push_back(static_cast<int>(i));
        const double cross_mean = replica_mean(r.cross_replicas);
        r.cross_sigma = std::sqrt(jackknife_entry(r.cross_replicas, cross_mean,
                                                  r.cross_replicas, cross_mean));
    }
    for (std::size_t row = 0; row < cov.fit_bins.size(); ++row) {
        const int i = cov.fit_bins[row];
        const auto &a = results[i].tangential_replicas;
        for (std::size_t col = diagonal_only ? row : 0; col <= row; ++col) {
            const int j = cov.fit_bins[col];
            const double v = jackknife_entry(a, means[i], results[j].tangential_replicas, means[j]);
            cov.matrix[i*n + j] = cov.matrix[j*n + i] = v;
        }
        results[i].sigma = std::sqrt(cov.matrix[i*n + i]);
    }
    return cov;
}

} // namespace bcg_esd
