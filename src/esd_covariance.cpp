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

// ==========================================
// Function: output_stream
// Method: Start a precise text output table and report an unusable destination immediately.
// ==========================================
std::ofstream output_stream(const std::string &path) {
    std::ofstream out(path, std::ios::trunc);
    if (!out) throw std::runtime_error("Cannot open output: " + path);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out << std::scientific << std::setprecision(16);
    return out;
}
} // namespace

// ==========================================
// Function: jackknife_covariance
// Method: Retain every regional deletion for all fitted radii and compute the full matrix.
// ==========================================
Covariance jackknife_covariance(std::vector<Result> &results, int regions) {
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
    for (const int i : cov.fit_bins) {
        const auto &a = results[i].tangential_replicas;
        for (const int j : cov.fit_bins) {
            if (j > i) break;
            const double v = jackknife_entry(a, means[i], results[j].tangential_replicas, means[j]);
            cov.matrix[i*n + j] = cov.matrix[j*n + i] = v;
        }
        results[i].sigma = std::sqrt(cov.matrix[i*n + i]);
    }
    return cov;
}

// ==========================================
// Function: write_results
// Method: Export measurements, ordered regional replicas, and the covariance used by the fit.
// ==========================================
void write_results(const std::vector<Result> &results, const Covariance &covariance,
                   const LensStack &stack, const Config &cfg) {
    auto out = output_stream(cfg.output);
    out << "# R_low R_high R_eff n_pair ESD_meas ESD_sigma ESD_cross ESD_cross_sigma\n"
        << "# Comoving R [Mpc/h], ESD [h M_sun/pc^2]\n"
        << "# center=" << bcg::center_name(cfg.center) << " groups=" << stack.n_lenses << " z_eff=" << stack.redshift
        << " mean_Mstar_Msun_over_h=" << stack.stellar_mass << " n_jack=" << stack.regions << '\n'
        << "# FD restored: " << cfg.fd_corrected << ", calibration: " << cfg.calibration.description() << '\n';
    if (cfg.self_calibrate)
        out << "# Self-calibration: Cartesian unique all-R union, sources=" << cfg.self_calibration_sources
            << ", FD interval=[" << -cfg.fd_limit() << ',' << cfg.fd_limit() << "], bins=" << cfg.selfcal_fd_bins << '\n';
    for (const auto &r : results)
        out << r.low << ' ' << r.high << ' ' << r.effective << ' ' << r.count << ' '
            << r.measured << ' ' << r.sigma << ' ' << r.cross << ' ' << r.cross_sigma << '\n';
    out.close();

    auto cov = output_stream(cfg.covariance_output);
    cov << "# Covariance [(h M_sun/pc^2)^2], rows and columns in ESD radial-bin order\n"
        << "# center=" << bcg::center_name(cfg.center) << " n_jack=" << stack.regions << " fit_bins_zero_based=";
    for (const int i : covariance.fit_bins) cov << ' ' << i;
    cov << '\n';
    const auto n = results.size();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            cov << covariance.matrix[i*n + j] << (j + 1 == n ? '\n' : ' ');
    cov.close();

    auto replicas = output_stream(cfg.replica_output);
    replicas << "# omitted_region ESD_tangential_R1..RN ESD_cross_R1..RN\n"
             << "# center=" << bcg::center_name(cfg.center) << '\n';
    for (int k = 0; k < stack.regions; ++k) {
        replicas << k;
        for (const auto &r : results) replicas << ' ' << r.tangential_replicas[k];
        for (const auto &r : results) replicas << ' ' << r.cross_replicas[k];
        replicas << '\n';
    }
    replicas.close();
}
} // namespace bcg_esd
