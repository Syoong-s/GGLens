#include "selfcalibration.hpp"
#include "bcg_archive.hpp"
#include <iostream>

namespace bcg_esd {
namespace {
void write_fd_points(const bcg::fs::path &path, const std::vector<double> &edges,
                     const std::array<std::vector<FdCalibrationPoint>, 2> &points,
                     std::int64_t selected) {
    bcg::archive::publish(path, [&](const auto &temporary) {
        std::ofstream out(temporary);
        if (!out) throw std::runtime_error("Cannot write FD points: " + path.string());
        out.exceptions(std::ios::badbit | std::ios::failbit);
        out << fd_points_header << '\n'
            << "# component 1 = g1/gf1, 2 = g2/gf2\n"
            << "# Bin numbers start at 1. FD intervals use gf_low <= gf < gf_high.\n"
            << "# status 0 = valid, 1 = too few sources, 2 = invalid PDF_SYM fit\n"
            << "# sigma is the PDF_SYM curvature error\n"
            << "# all_R_sources=" << selected << '\n'
            << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10 - 1);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t component = 0; component < points.size(); ++component) {
            for (std::size_t bin = 0; bin < points[component].size(); ++bin) {
                const auto &p = points[component][bin];
                const bool valid = p.status == 0;
                out << component + 1 << ' ' << bin + 1 << ' '
                    << edges[bin] << ' ' << edges[bin + 1] << ' ' << (p.count ? p.field : nan) << ' '
                    << (valid ? p.value : nan) << ' ' << (valid ? p.sigma : nan) << ' '
                    << p.count << ' ' << p.status << '\n';
            }
        }
        out.close();
    });
}
} // namespace

// ==========================================
// Function: fit_selfcal_line
// Method: Fit g=(1+m)*gf+c with absolute-sigma weights and centered, scaled sums.
// ==========================================
FdLineFit fit_selfcal_line(const std::vector<FdCalibrationPoint> &points) {
    double min_sigma = std::numeric_limits<double>::infinity();
    for (const auto &p : points) if (p.status == 0 && p.count > 0 &&
        std::isfinite(p.field) && std::isfinite(p.value) && std::isfinite(p.sigma) && p.sigma > 0)
        min_sigma = std::min(min_sigma, p.sigma);
    if (!std::isfinite(min_sigma)) throw std::runtime_error("Self-calibration has no valid FD points");
    long double weight_sum = 0, x_sum = 0, y_sum = 0;
    FdLineFit fit;
    for (const auto &p : points) if (p.status == 0 && p.count > 0 &&
        std::isfinite(p.field) && std::isfinite(p.value) && std::isfinite(p.sigma) && p.sigma > 0) {
        const long double w = std::pow(static_cast<long double>(min_sigma) / p.sigma, 2);
        weight_sum += w; x_sum += w*p.field; y_sum += w*p.value; ++fit.count;
    }
    if (fit.count < 2) throw std::runtime_error("Self-calibration needs at least two valid FD bins per component");
    const long double x_mean = x_sum/weight_sum, y_mean = y_sum/weight_sum;
    long double spread = 0, product = 0;
    for (const auto &p : points) if (p.status == 0 && p.count > 0 &&
        std::isfinite(p.field) && std::isfinite(p.value) && std::isfinite(p.sigma) && p.sigma > 0) {
        const long double w = std::pow(static_cast<long double>(min_sigma) / p.sigma, 2);
        const long double dx = p.field-x_mean;
        spread += w*dx*dx; product += w*dx*(p.value-y_mean);
    }
    if (!(spread > 0) || !std::isfinite(spread))
        throw std::runtime_error("Self-calibration FD bins have no finite field-distortion spread");
    const long double slope = product/spread;
    fit.m = static_cast<double>(slope-1); fit.c = static_cast<double>(y_mean-slope*x_mean);
    fit.sigma_m = static_cast<double>(min_sigma/std::sqrt(spread));
    fit.sigma_c = static_cast<double>(min_sigma*std::sqrt(1/weight_sum+x_mean*x_mean/spread));
    fit.covariance = static_cast<double>(-min_sigma*min_sigma*x_mean/spread);
    for (const auto &p : points) if (p.status == 0 && p.count > 0 &&
        std::isfinite(p.field) && std::isfinite(p.value) && std::isfinite(p.sigma) && p.sigma > 0)
        fit.chi2 += std::pow((p.value-(1+fit.m)*p.field-fit.c)/p.sigma, 2);
    if (!(std::isfinite(fit.m) && std::isfinite(fit.c) && std::isfinite(fit.sigma_m) && std::isfinite(fit.sigma_c)))
        throw std::runtime_error("Nonfinite self-calibration line fit");
    return fit;
}

// ==========================================
// Function: SelfCalibration constructor
// Method: Build equal-width FD bins and allocate their spools within the rank's memory quota.
// ==========================================
SelfCalibration::SelfCalibration(const Config &cfg, const bcg::ScratchSpace &scratch)
    : cfg_(cfg), scratch_(scratch), edges_(static_cast<std::size_t>(cfg.selfcal_fd_bins)+1) {
    const double limit = cfg.fd_limit();
    const double width = 2 * limit / cfg.selfcal_fd_bins;
    for (int i = 0; i <= cfg.selfcal_fd_bins; ++i) edges_[i] = -limit + width * i;
    edges_.front() = -limit;
    edges_.back() = limit;
    const auto bytes = cfg.memory_bytes()/16/(2*static_cast<std::size_t>(cfg.selfcal_fd_bins));
    for (int component = 0; component < 2; ++component) {
        field_sum_[component].assign(cfg.selfcal_fd_bins, 0);
        for (int bin = 0; bin < cfg.selfcal_fd_bins; ++bin)
            cells_[component].emplace_back(new bcg::PairStore(
                scratch.path("selfcal_" + std::to_string(component) + "_" + std::to_string(bin) + ".bin"),
                1024*1024, bytes, 1, 1, true));
    }
}

// ==========================================
// Function: SelfCalibration::append
// Method: Accept each matched raw measurement row once, without lens multiplicity or rotation.
// ==========================================
void SelfCalibration::append(const bcg::Source &source) {
    if (local_selected_ == std::numeric_limits<std::int64_t>::max())
        throw std::runtime_error("Self-calibration source count exceeds int64");
    ++local_selected_;
    for (int component = 0; component < 2; ++component) {
        const double field = component ? source.gf2 : source.gf1;
        if (field < edges_.front() || field >= edges_.back()) continue;
        const auto bin = static_cast<std::size_t>(std::upper_bound(edges_.begin(), edges_.end(), field)-edges_.begin()-1);
        const double g = component ? source.g2 : source.g1;
        const double response = component ? source.de+source.h1 : source.de-source.h1;
        cells_[component][bin]->append({g, 0, response, 0, 0, 0}, 1, 0);
        field_sum_[component][bin] += field;
    }
}

// ==========================================
// Function: SelfCalibration::finish
// Method: Save the Cartesian all-R bin measurements before fitting and validating both lines.
// ==========================================
bcg::ShearCalibration SelfCalibration::finish(MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    MPI_Allreduce(&local_selected_, &total_selected_, 1, MPI_INT64_T, MPI_SUM, comm);
    if (!total_selected_) throw std::runtime_error("Self-calibration has no selected all-R source rows");
    std::array<std::vector<FdCalibrationPoint>, 2> points;
    for (auto &component : cells_) for (auto &cell : component) cell->finish();
    for (int component = 0; component < 2; ++component) {
        points[component].resize(cfg_.selfcal_fd_bins);
        for (int bin = 0; bin < cfg_.selfcal_fd_bins; ++bin) {
            auto &cell = *cells_[component][bin]; auto &point = points[component][bin];
            MPI_Allreduce(&cell.count, &point.count, 1, MPI_INT64_T, MPI_SUM, comm);
            double sum = 0;
            MPI_Allreduce(&field_sum_[component][bin], &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
            if (point.count) point.field = sum/point.count;
            if (point.count >= cfg_.selfcal_min_count) {
                cell.load_if_small(cfg_.memory_bytes()/4);
                const auto fit = selfcal_fd_pdf(cell, cfg_.selfcal_pdf_bins, cfg_.selfcal_min_count,
                    cfg_.selfcal_search_width, scratch_, cfg_.memory_bytes()/8, comm);
                point.status = fit.valid ? 0 : 2;
                if (fit.valid) { point.value = fit.value; point.sigma = fit.sigma; }
            }
            if (!rank) std::cout << "Self-calibration g" << component+1 << " FD bin " << bin+1 << '/'
                << cfg_.selfcal_fd_bins << " rows=" << point.count << " status=" << point.status << std::endl;
            cell.release();
        }
    }
    std::array<double, 4> coefficients{};
    if (!rank) {
        write_fd_points(cfg_.fd_points_output, edges_, points, total_selected_);
        std::cout << "Saved FD calibration points: " << cfg_.fd_points_output << std::endl;
        const std::array<FdLineFit, 2> fits{fit_selfcal_line(points[0]), fit_selfcal_line(points[1])};
        bcg::ShearCalibration fitted{fits[0].m, fits[1].m, fits[0].c, fits[1].c};
        fitted.validate();
        coefficients = {fitted.m1, fitted.m2, fitted.c1, fitted.c2};
        std::cout << "Self-calibration all-R unique rows=" << total_selected_ << "; "
                  << fitted.description() << std::endl;
    }
    MPI_Bcast(coefficients.data(), 4, MPI_DOUBLE, 0, comm);
    bcg::ShearCalibration result{coefficients[0], coefficients[1], coefficients[2], coefficients[3]};
    result.validate(); return result;
}
} // namespace bcg_esd
