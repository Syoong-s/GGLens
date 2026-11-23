#pragma once
#include "bcg_catalog.hpp"

namespace bcg {
struct Cuts;

// ==========================================
// Structure: ExternalColumns
// Method: Store the configured positions of catalog metadata.
// ==========================================
struct ExternalColumns {
    int ra, dec, zp, mag, ext_mash, expo_num, ccd_num, expo_chi2;
};

// ==========================================
// Structure: ShearRow
// Method: Give each consumed catalog field a stable slot for MPI column mapping.
// ==========================================
struct ShearRow {
    enum Field {
        ra, dec, zp, mag, ext_mash, expo_num, ccd_num, chi2,
        g1, g2, de, h1, h2, gf1, gf2,
        snr_f, half_light_flux, half_light_area, psf_fwhm,
        nstar, poly_chi2, flag, sigma, imax, jmax, count
    };
    std::array<double, count> values{};
};

// ==========================================
// Structure: ShearLayout
// Method: Store the first catalog width and the destination of each input column.
// ==========================================
struct ShearLayout {
    std::vector<int> destinations;
};

ShearLayout read_shear_layout(const std::string &first_path, const ExternalColumns &external,
                              const Cuts &cuts, MPI_Comm comm);
void for_each_shear_row(const std::string &path, const ShearLayout &layout,
                        const std::function<void(const ShearRow &)> &visit);
} // namespace bcg
