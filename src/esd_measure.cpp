#include "bcg_esd.hpp"
#include "bcg_physics.hpp"
#include "bcg_pdf_stream.hpp"
#include "bcg_quantiles.hpp"
#include <iostream>

namespace bcg_esd {
// ==========================================
// Function: measure_bin
// Method: Measure bounded pair views with NFW-centered initial search and jackknife statistics.
// ==========================================
Result measure_bin(bcg::PairStore &pairs, double low, double high,
                   const std::vector<bcg::Lens> &lenses, const std::vector<int> &labels,
                   const Config &cfg, const bcg::ScratchSpace &scratch, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Result out{low, high, nan, 0, nan, nan, nan, nan, nan, nan, nan};
    MPI_Allreduce(&pairs.count, &out.count, 1, MPI_INT64_T, MPI_SUM, comm);
    if (out.count < static_cast<std::int64_t>(cfg.pdf_bins) + 2) return out;
    const bool resident = pairs.load_if_small(cfg.memory_bytes() / 2);
    int local_streaming = !resident && pairs.count > 0 ? 1 : 0, streaming = 0;
    MPI_Allreduce(&local_streaming, &streaming, 1, MPI_INT, MPI_SUM, comm);
    double global_logr = 0;
    MPI_Allreduce(&pairs.log_radius_sum, &global_logr, 1, MPI_DOUBLE, MPI_SUM, comm);
    std::vector<unsigned char> present(lenses.size(), 0);
    MPI_Allreduce(pairs.present.data(), present.data(), static_cast<int>(present.size()),
                  MPI_UNSIGNED_CHAR, MPI_MAX, comm);
    out.effective = std::exp(global_logr / out.count);
    double theory = 0;
    if (!rank) {
        int active = 0;
        for (std::size_t i = 0; i < lenses.size(); ++i) if (present[i]) {
            const double mass = std::pow(10.0, lenses[i].log_mh);
            const double conc = bcg::concentration(mass, lenses[i].z, cfg.conc_model,
                                                  cfg.nfw_c, cfg.conc_min, cfg.conc_max);
            theory += bcg::nfw_delta_sigma(out.effective, mass, conc, cfg.omega_m); ++active;
        }
        if (!active) throw std::runtime_error("Pairs exist but no active lens found");
        theory /= active;
    }
    MPI_Bcast(&theory, 1, MPI_DOUBLE, 0, comm); out.theory = theory;
    const double width = std::max(5.0 * std::abs(theory), 10.0);
    std::vector<std::int64_t> omitted(cfg.n_jack, 0);
    MPI_Allreduce(pairs.region_count.data(), omitted.data(), cfg.n_jack, MPI_INT64_T, MPI_SUM, comm);
    std::vector<double> weights(cfg.n_jack, 0);
    for (int r = 0; r < cfg.n_jack; ++r)
        if (omitted[r] > 0 && out.count - omitted[r] >= static_cast<std::int64_t>(cfg.pdf_bins) + 2)
            weights[r] = static_cast<double>(omitted[r]);
    if (!rank) std::cout << "  streaming workers=" << streaming << ", pairs=" << out.count << std::endl;
    for (int component = 0; component < 2; ++component) {
        if (!rank) std::cout << "  " << (component ? "cross" : "tangential")
                             << ": exact quantiles and full PDF" << std::endl;
        const auto initial = bcg::prepare_quantiles(pairs, labels, component, cfg.pdf_bins,
                                                    cfg.n_jack, scratch, cfg.memory_bytes() / 8);
        const double center = component ? 0.0 : theory;
        const bcg::PairView full_view{pairs, labels, component, -1};
        const auto full = bcg::pdf_sym_stream(full_view, initial.values[0], cfg.pdf_bins, comm, center, width);
        std::vector<double> estimates(cfg.n_jack, 0);
        for (int r = 0; r < cfg.n_jack; ++r) {
            if (weights[r] <= 0) continue;
            const bcg::PairView view{pairs, labels, component, r};
            estimates[r] = bcg::pdf_sym_stream(view, initial.values[static_cast<std::size_t>(r) + 1],
                cfg.pdf_bins, comm, center, width, out.count).value;
            if (!rank && ((r + 1) % 10 == 0 || r + 1 == cfg.n_jack))
                std::cout << "    jackknife " << r + 1 << '/' << cfg.n_jack << std::endl;
        }
        const double sigma = bcg::jackknife_sigma(estimates, weights, full.value, out.count);
        const double error = std::isfinite(sigma) && sigma > 0 ? sigma : full.fit_sigma;
        if (!component) { out.measured = full.value; out.sigma = error; }
        else { out.cross = full.value; out.cross_sigma = error; }
    }
    out.ratio = theory != 0 ? out.measured / theory : nan;
    out.residual = out.measured - theory; return out;
}

// ==========================================
// Function: write_results
// Method: Export eleven ESD columns with effective manual/automatic coefficients and all-R FD provenance.
// ==========================================
void write_results(const std::vector<Result> &results, const Config &cfg) {
    const auto &path = cfg.output;
    std::ofstream out(path, std::ios::trunc);
    if (!out) throw std::runtime_error("Cannot open output: " + path);
    out << "# R_low R_high R_eff n_pair ESD_meas ESD_sigma ESD_cross ESD_cross_sigma ESD_true ratio residual\n"
        << "# Units: R in h^-1 Mpc, ESD in h M_sun/pc^2\n"
        << "# BCG-centered PDF_SYM; M_180m NFW reference\n"
        << "# Catalog FD restored: " << cfg.fd_corrected << '\n'
        << "# Shear calibration: " << cfg.calibration.description() << '\n'
        << "# G1_corr=(G1-c1*(de-h1))/(1+m1)-gf1*(de-h1)+gf2*h2\n"
        << "# G2_corr=(G2-c2*(de+h1))/(1+m2)-gf2*(de+h1)+gf1*h2\n"
        << "# ESD errors condition on fixed m1/m2/c1/c2; calibration uncertainty not included\n";
    if (cfg.self_calibrate)
        out << "# Self-calibration: Cartesian unique all-R union; sources=" << cfg.self_calibration_sources
            << " fd_bins=" << cfg.selfcal_fd_bins << " magnitude_cells=" << cfg.selfcal_pdf_bins
            << " plot=" << cfg.fd_plot << '\n';
    out << std::scientific << std::setprecision(8);
    for (const auto &r : results)
        out << r.low << ' ' << r.high << ' ' << r.effective << ' ' << r.count << ' '
            << r.measured << ' ' << r.sigma << ' ' << r.cross << ' ' << r.cross_sigma << ' '
            << r.theory << ' ' << r.ratio << ' ' << r.residual << '\n';
    out.flush();
    if (!out) throw std::runtime_error("Failed writing output: " + path);
}

} // namespace bcg_esd
