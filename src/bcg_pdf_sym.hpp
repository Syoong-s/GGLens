#pragma once

#include <mpi.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace bcg {
struct PdfResult { double value = 0, fit_sigma = 0, chi2 = 0; };

// ==========================================
// Function: f77_bin
// Method: First upper edge >= |x| wins, including nonmonotonic refinements.
// ==========================================
inline int f77_bin(double value, const std::vector<double> &edges) {
    for (int k = 0; k < static_cast<int>(edges.size()); ++k)
        if (value <= edges[k]) return k;
    return static_cast<int>(edges.size());
}

// ==========================================
// Function: solve_quadratic
// Method: Scaled least-squares quadratic fit through a 3x3 normal system.
// ==========================================
inline std::array<double, 3> solve_quadratic(const std::vector<double> &c,
                                             const std::vector<double> &y) {
    if (c.size() < 3) return {0, 0, 0};
    const double center = c[c.size()/2];
    double scale = c.back() - c.front();
    if (scale == 0) scale = 1;
    double s[5] = {}, t[3] = {};
    for (std::size_t i = 0; i < c.size(); ++i) {
        const double u = (c[i] - center) / scale;
        const double u2 = u*u;
        s[0] += 1; s[1] += u; s[2] += u2; s[3] += u2*u; s[4] += u2*u2;
        t[0] += y[i]; t[1] += u*y[i]; t[2] += u2*y[i];
    }
    double a[3][4] = {{s[4], s[3], s[2], t[2]},
                      {s[3], s[2], s[1], t[1]},
                      {s[2], s[1], s[0], t[0]}};
    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        for (int row = col+1; row < 3; ++row)
            if (std::abs(a[row][col]) > std::abs(a[pivot][col])) pivot = row;
        if (std::abs(a[pivot][col]) < 1e-20) return {0, 0, 0};
        for (int j = col; j < 4; ++j) std::swap(a[col][j], a[pivot][j]);
        const double div = a[col][col];
        for (int j = col; j < 4; ++j) a[col][j] /= div;
        for (int row = 0; row < 3; ++row) if (row != col) {
            const double f = a[row][col];
            for (int j = col; j < 4; ++j) a[row][j] -= f * a[col][j];
        }
    }
    const double a1 = a[0][3] / (scale*scale);
    const double a2 = a[1][3] / scale - 2*a1*center;
    const double a3 = a[2][3] - a1*center*center - a2*center;
    return {a1, a2, a3};
}

// ==========================================
// Function: jackknife_sigma
// Method: Reproduce the weighted leave-one-region-out variance expression.
// ==========================================
inline double jackknife_sigma(const std::vector<double> &estimates,
                              const std::vector<double> &weights,
                              double full, double ntot) {
    int used = 0;
    double exterior = 0;
    for (std::size_t i = 0; i < weights.size(); ++i) if (weights[i] > 0) {
        ++used; exterior += (1.0 - weights[i]/ntot) * estimates[i];
    }
    if (used < 2) return 0;
    double sum = 0;
    for (std::size_t i = 0; i < weights.size(); ++i) if (weights[i] > 0) {
        const double h = ntot / weights[i];
        const double t = (h - used)*full - (h - 1)*estimates[i] + exterior;
        sum += t*t / (h - 1);
    }
    return std::sqrt(sum / used);
}
} // namespace bcg
