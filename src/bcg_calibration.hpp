#pragma once
#include "bcg_catalog.hpp"

namespace bcg {
// ==========================================
// Function: restore_raw_shear
// Method: Invert the catalog's Cartesian FD subtraction, including both h2 cross terms.
// ==========================================
inline void restore_raw_shear(Source &source) {
    source.g1 += source.gf1 * (source.de - source.h1) - source.gf2 * source.h2;
    source.g2 += source.gf2 * (source.de + source.h1) - source.gf1 * source.h2;
    if (!(std::isfinite(source.g1) && std::isfinite(source.g2)))
        throw std::runtime_error("Non-finite raw G1/G2 after restoring catalog FD");
}

// ==========================================
// Structure: ShearCalibration
// Method: Store independent Cartesian FD-fit biases; field distortion is always removed.
// ==========================================
struct ShearCalibration {
    double m1 = 0, m2 = 0, c1 = 0, c2 = 0;

    // ==========================================
    // Function: validate
    // Method: Require finite FD-fit coefficients and positive multiplicative gains.
    // ==========================================
    void validate() const {
        if (!(std::isfinite(m1) && std::isfinite(m2) &&
              std::isfinite(c1) && std::isfinite(c2) && m1 > -1.0 && m2 > -1.0))
            throw std::runtime_error("Calibration requires finite m1/m2/c1/c2 and m1,m2 > -1");
    }

    // ==========================================
    // Function: apply
    // Method: Subtract independent FD-fit intercepts, divide gains, then remove the full FD response.
    // ==========================================
    void apply(Source &source) const {
        const double b1 = source.de - source.h1;
        const double b2 = source.de + source.h1;
        const double g1 = (source.g1 - c1 * b1) / (1.0 + m1)
                        - source.gf1 * b1 + source.gf2 * source.h2;
        const double g2 = (source.g2 - c2 * b2) / (1.0 + m2)
                        - source.gf2 * b2 + source.gf1 * source.h2;
        if (!(std::isfinite(g1) && std::isfinite(g2)))
            throw std::runtime_error("Non-finite calibrated G1/G2; check calibration coefficients and catalog moments");
        source.g1 = g1;
        source.g2 = g2;
    }

    // ==========================================
    // Function: description
    // Method: Serialize calibration coefficients and conventions for logs and result provenance.
    // ==========================================
    std::string description() const {
        std::ostringstream out;
        out << std::setprecision(17) << "m1=" << m1 << " m2=" << m2
            << " c1=" << c1 << " c2=" << c2
            << " frame=celestial c_model=component_fd fd_subtracted=1";
        return out.str();
    }
};
} // namespace bcg
