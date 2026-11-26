#pragma once
#include "bcg_shear_catalog.hpp"
#include <optional>

namespace bcg {
// ==========================================
// Function: cut_enabled
// Method: Honor the disabled threshold for scalar catalog cuts.
// ==========================================
constexpr bool cut_enabled(double threshold) {
    return threshold >= -999.0;
}

// ==========================================
// Structure: CutThresholds
// Method: Store the scalar cuts for compile-time initialization and INI overrides.
// ==========================================
struct CutThresholds {
    double snrfcut, snrlow, snrhigh, starcut;
    double chi2_thresh, flagcut, zplow, zphigh;
    double gf_lim, imaxcut, jmaxcut;
    double psf_chi2_mltp, r_half_thresh;
    double star_bar_mltp, expo_bar_mltp, max_dup;
};

struct ExternalColumnCut {
    std::string key;
    int column;
    std::optional<double> low, high;
};

// ==========================================
// Structure: Cuts
// Method: Keep the configured cuts and their shared column dependencies together.
// ==========================================
struct Cuts : CutThresholds {
    std::vector<int> excluded_ccds; // Sorted during configuration validation.
    std::vector<ExternalColumnCut> external{};

    // ==========================================
    // Function: ccd_cut_enabled
    // Method: Require the CCD column only when at least one ID is excluded.
    // ==========================================
    bool ccd_cut_enabled() const { return !excluded_ccds.empty(); }

    // ==========================================
    // Function: exposure_bar_enabled / bar_enabled
    // Method: Use the common cut switch and keep per-exposure precedence.
    // ==========================================
    bool exposure_bar_enabled() const { return cut_enabled(expo_bar_mltp); }
    bool bar_enabled() const { return cut_enabled(star_bar_mltp) || exposure_bar_enabled(); }

    // ==========================================
    // Function: needs_snr / needs_psf
    // Method: Identify the columns consumed by the enabled cuts and stellar-bar fit.
    // ==========================================
    bool needs_snr() const { return bar_enabled() || cut_enabled(snrlow) || cut_enabled(snrhigh); }
    bool needs_psf() const { return bar_enabled() || cut_enabled(r_half_thresh); }
};

void for_each_local_source(const std::vector<std::string> &paths, const Cuts &cuts,
                           const ShearLayout &layout, MPI_Comm comm,
                           const std::function<void(const Source &, const BarDatum &)> &visit);
LocalSources read_local_sources(const std::vector<std::string> &paths, const Cuts &cuts,
                                const ShearLayout &layout, MPI_Comm comm);
} // namespace bcg
