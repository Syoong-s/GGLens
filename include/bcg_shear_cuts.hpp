#pragma once
#include "bcg_shear_catalog.hpp"

namespace bcg {
// ==========================================
// Function: cut_enabled
// Method: Honor the common disabled threshold for all catalog cuts.
// ==========================================
constexpr bool cut_enabled(double threshold) {
    return threshold >= -999.0;
}

// ==========================================
// Structure: Cuts
// Method: Keep input thresholds and their shared column dependencies together.
// ==========================================
struct Cuts {
    double snrfcut, snrlow, snrhigh, starcut;
    double chi2_thresh, flagcut, zplow, zphigh;
    double gf_lim, imaxcut, jmaxcut;
    double mag_low, mag_high, psf_chi2_mltp;
    double r_half_thresh, ext_cut;
    double star_bar_mltp, expo_bar_mltp, max_dup;

    // ==========================================
    // Function: exposure_bar_enabled / bar_enabled
    // Method: Use the common cut switch and keep per-exposure precedence.
    // ==========================================
    bool exposure_bar_enabled() const { return cut_enabled(expo_bar_mltp); }
    bool bar_enabled() const { return cut_enabled(star_bar_mltp) || exposure_bar_enabled(); }

    // ==========================================
    // Function: needs_magnitude / needs_snr / needs_psf
    // Method: Identify the columns consumed by the enabled cuts and stellar-bar fit.
    // ==========================================
    bool needs_magnitude() const { return bar_enabled() || cut_enabled(mag_low) || cut_enabled(mag_high); }
    bool needs_snr() const { return bar_enabled() || cut_enabled(snrlow) || cut_enabled(snrhigh); }
    bool needs_psf() const { return bar_enabled() || cut_enabled(r_half_thresh); }
};

void for_each_local_source(const std::vector<std::string> &paths, const Cuts &cuts,
                           const ShearLayout &layout, MPI_Comm comm,
                           const std::function<void(const Source &, const BarDatum &)> &visit);
LocalSources read_local_sources(const std::vector<std::string> &paths, const Cuts &cuts,
                                const ShearLayout &layout, MPI_Comm comm);
} // namespace bcg
