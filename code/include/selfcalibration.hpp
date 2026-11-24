#pragma once
#include "bcg_esd.hpp"

namespace bcg_esd {
inline constexpr char fd_points_header[] = "# component bin gf_low gf_high gf_mean shear sigma n_sources status";

// ==========================================
// Structure: FdCalibrationPoint
// Method: Retain one global Cartesian FD-bin estimate and its fit validity.
// ==========================================
struct FdCalibrationPoint {
    std::int64_t count = 0;
    double field = 0, value = 0, sigma = 0;
    int status = 1;
};

// ==========================================
// Structure: FdLineFit
// Method: Describe the absolute-sigma weighted g=(1+m)*gf+c line.
// ==========================================
struct FdLineFit {
    double m = 0, c = 0, sigma_m = 0, sigma_c = 0, covariance = 0, chi2 = 0;
    int count = 0;
};

// ==========================================
// Structure: FdPdfFit
// Method: Return the original cochko zero-centered PDF fit or an invalid point.
// ==========================================
struct FdPdfFit {
    double value = 0, sigma = 0;
    bool valid = false;
};
FdPdfFit selfcal_fd_pdf(const bcg::PairStore &store, int magnitude_bins, int min_count,
                       double search_width, const bcg::ScratchSpace &scratch,
                       std::size_t sort_bytes, MPI_Comm comm);
FdLineFit fit_selfcal_line(const std::vector<FdCalibrationPoint> &points);

// ==========================================
// Class: SelfCalibration
// Method: Spool unique raw all-R rows and fit the two Cartesian FD components.
// ==========================================
class SelfCalibration {
    const Config &cfg_;
    const bcg::ScratchSpace &scratch_;
    std::vector<double> edges_;
    std::array<std::vector<std::unique_ptr<bcg::PairStore>>, 2> cells_;
    std::array<std::vector<double>, 2> field_sum_;
    std::int64_t local_selected_ = 0, total_selected_ = 0;
public:
    SelfCalibration(const Config &cfg, const bcg::ScratchSpace &scratch);
    void append(const bcg::Source &source);
    bcg::ShearCalibration finish(MPI_Comm comm);
    // ==========================================
    // Function: selected_sources
    // Method: Expose the globally reduced unique all-R row count for ESD provenance.
    // ==========================================
    std::int64_t selected_sources() const { return total_selected_; }
};
} // namespace bcg_esd
