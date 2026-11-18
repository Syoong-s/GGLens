#pragma once

#include "bcg_catalog.hpp"
#include <array>
#include <numeric>
#include <random>
#include <utility>

namespace bcg {
constexpr double pi = 3.1415926535897932384626433832795;
constexpr double deg = pi / 180.0;

// ==========================================
// Function: sky_xyz
// Method: Convert equatorial coordinates to a unit-vector for spatial search.
// ==========================================
inline std::array<double, 3> sky_xyz(double ra, double dec) {
    const double a = ra * deg, d = dec * deg, c = std::cos(d);
    return {c * std::cos(a), c * std::sin(a), std::sin(d)};
}

class Cosmology {
    double omega_m_, h_, dz_;
    std::vector<double> chi_;
public:
    // ==========================================
    // Function: Cosmology constructor
    // Method: Precompute comoving distance using Simpson integration per cell.
    // ==========================================
    Cosmology(double omega_m, double h, double zmax = 5.0, int n = 4001)
        : omega_m_(omega_m), h_(h), dz_(zmax / (n - 1)), chi_(n, 0.0) {
        for (int i = 1; i < n; ++i) {
            const double a = (i - 1) * dz_, b = i * dz_;
            chi_[i] = chi_[i - 1] + 2998.0 * dz_ *
                (inv_E(a) + 4.0 * inv_E((a + b) * 0.5) + inv_E(b)) / 6.0;
        }
    }

    // ==========================================
    // Function: inv_E
    // Method: Flat Lambda-CDM 1/E(z).
    // ==========================================
    double inv_E(double z) const {
        return 1.0 / std::sqrt(omega_m_ * std::pow(1.0 + z, 3) + 1.0 - omega_m_);
    }

    // ==========================================
    // Function: chi
    // Method: Cubic Hermite interpolation with exact integrand end slopes.
    // ==========================================
    double chi(double z) const {
        if (z <= 0) return 0;
        const double maxz = dz_ * (chi_.size() - 1);
        if (z > maxz) throw std::runtime_error("Redshift exceeds cosmology table");
        const int i = std::min(static_cast<int>(z / dz_), static_cast<int>(chi_.size()) - 2);
        const double t = (z - i * dz_) / dz_;
        const double y0 = chi_[i], y1 = chi_[i + 1];
        const double m0 = 2998.0 * inv_E(i * dz_);
        const double m1 = 2998.0 * inv_E((i + 1) * dz_);
        return (2*t*t*t - 3*t*t + 1)*y0 + (t*t*t - 2*t*t + t)*dz_*m0 +
               (-2*t*t*t + 3*t*t)*y1 + (t*t*t - t*t)*dz_*m1;
    }

    // ==========================================
    // Function: beta
    // Method: Inverse comoving Sigma_crit in pc^2/(h M_sun).
    // ==========================================
    double beta(double zl, double zs, double chil, double chis) const {
        if (zs <= zl || chis <= chil || chil <= 0) return 0;
        return 1e12 * chil * (chis - chil) * (1.0 + zl) / (1.661e18 * chis);
    }
};

class KDTree {
    struct Node {
        int point = -1, left = -1, right = -1, axis = 0;
        std::array<double, 3> lo{}, hi{};
    };
    const std::vector<Source> &sources_;
    std::vector<int> order_;
    std::vector<Node> nodes_;
    int root_ = -1;

    // ==========================================
    // Function: KDTree::build
    // Method: Median split the three-dimensional unit-sphere coordinates.
    // ==========================================
    int build(int begin, int end, int depth) {
        if (begin >= end) return -1;
        const int axis = depth % 3, mid = begin + (end - begin) / 2;
        std::nth_element(order_.begin() + begin, order_.begin() + mid,
                         order_.begin() + end, [&](int a, int b) {
            return sources_[a].xyz[axis] < sources_[b].xyz[axis];
        });
        const int node = static_cast<int>(nodes_.size());
        nodes_.push_back({});
        nodes_[node].point = order_[mid];
        nodes_[node].axis = axis;
        const int left = build(begin, mid, depth + 1);
        const int right = build(mid + 1, end, depth + 1);
        nodes_[node].left = left;
        nodes_[node].right = right;
        for (int k = 0; k < 3; ++k) {
            double lo = sources_[order_[mid]].xyz[k], hi = lo;
            if (left >= 0) { lo = std::min(lo, nodes_[left].lo[k]); hi = std::max(hi, nodes_[left].hi[k]); }
            if (right >= 0) { lo = std::min(lo, nodes_[right].lo[k]); hi = std::max(hi, nodes_[right].hi[k]); }
            nodes_[node].lo[k] = lo; nodes_[node].hi[k] = hi;
        }
        return node;
    }

    // ==========================================
    // Function: KDTree::query_node
    // Method: Prune by bounding box, then test exact chord distance.
    // ==========================================
    void query_node(int node, const std::array<double, 3> &q, double r2,
                    std::vector<int> &found) const {
        if (node < 0) return;
        const auto &n = nodes_[node];
        double box_d2 = 0, point_d2 = 0;
        for (int k = 0; k < 3; ++k) {
            if (q[k] < n.lo[k]) box_d2 += (n.lo[k] - q[k]) * (n.lo[k] - q[k]);
            else if (q[k] > n.hi[k]) box_d2 += (q[k] - n.hi[k]) * (q[k] - n.hi[k]);
            const double d = q[k] - sources_[n.point].xyz[k];
            point_d2 += d * d;
        }
        if (box_d2 > r2) return;
        if (point_d2 <= r2) found.push_back(n.point);
        query_node(n.left, q, r2, found);
        query_node(n.right, q, r2, found);
    }
public:
    // ==========================================
    // Function: KDTree constructor
    // Method: Index all local sources once for all selected lenses.
    // ==========================================
    explicit KDTree(const std::vector<Source> &sources) : sources_(sources) {
        order_.resize(sources.size());
        std::iota(order_.begin(), order_.end(), 0);
        nodes_.reserve(sources.size());
        root_ = build(0, static_cast<int>(sources.size()), 0);
    }

    // ==========================================
    // Function: KDTree::radius
    // Method: Return source indices within a spherical chord radius.
    // ==========================================
    void radius(const std::array<double, 3> &q, double radius,
                std::vector<int> &found) const {
        found.clear();
        query_node(root_, q, radius * radius, found);
    }
};

// ==========================================
// Function: concentration
// Method: Evaluate the requested c(M,z) law with finite clipping.
// ==========================================
inline double concentration(double mass, double z, const std::string &model,
                            double fixed, double cmin, double cmax) {
    double c = fixed;
    if (model == "duffy08_200m") c = 10.14 * std::pow(mass / 2e12, -0.081) * std::pow(1 + z, -1.01);
    else if (model == "duffy08_200c") c = 5.71 * std::pow(mass / 2e12, -0.084) * std::pow(1 + z, -0.47);
    else if (model == "duffy08_vir") c = 7.85 * std::pow(mass / 2e12, -0.081) * std::pow(1 + z, -0.71);
    else if (model == "dutton14_200c") {
        const double a = 0.520 + 0.385 * std::exp(-0.617 * std::pow(z, 1.21));
        const double b = -0.101 + 0.026 * z;
        c = std::pow(10.0, a + b * std::log10(mass / 1e12));
    } else if (model != "fixed") throw std::runtime_error("Unknown concentration model: " + model);
    if (!std::isfinite(c)) c = fixed;
    return std::clamp(c, cmin, cmax);
}

// ==========================================
// Function: nfw_f
// Method: Wright-Brainerd projected NFW surface-density kernel.
// ==========================================
inline double nfw_f(double x) {
    if (std::abs(x - 1.0) <= 1e-10) return 1.0 / 3.0;
    if (x < 1.0) {
        const double s = std::sqrt(1.0 - x*x);
        return (1.0 - std::acosh(1.0 / x) / s) / (x*x - 1.0);
    }
    const double s = std::sqrt(x*x - 1.0);
    return (1.0 - std::acos(1.0 / x) / s) / (x*x - 1.0);
}

// ==========================================
// Function: gauss_legendre_integral
// Method: Fixed 96-point quadrature; split at x=1 in the caller.
// ==========================================
inline double gauss_legendre_integral(double a, double b) {
    if (b <= a) return 0;
    static const auto table = [] {
        std::array<std::pair<double, double>, 96> q{};
        for (int i = 0; i < 48; ++i) {
            double z = std::cos(pi * (i + 0.75) / 96.5), prev;
            double pp = 0;
            do {
                double p1 = 1, p2 = 0;
                for (int j = 1; j <= 96; ++j) {
                    const double p3 = p2;
                    p2 = p1;
                    p1 = ((2*j - 1)*z*p2 - (j - 1)*p3) / j;
                }
                pp = 96.0 * (z*p1 - p2) / (z*z - 1);
                prev = z;
                z = prev - p1 / pp;
            } while (std::abs(z - prev) > 1e-15);
            const double w = 2.0 / ((1.0 - z*z)*pp*pp);
            q[i] = {-z, w}; q[95-i] = {z, w};
        }
        return q;
    }();
    const double half = (b - a) * 0.5, center = (a + b) * 0.5;
    double sum = 0;
    for (const auto &nw : table) {
        const double x = center + half * nw.first;
        sum += nw.second * x * nfw_f(x);
    }
    return half * sum;
}

// ==========================================
// Function: nfw_delta_sigma
// Method: M_180m NFW excess density in h M_sun/pc^2.
// ==========================================
inline double nfw_delta_sigma(double radius, double mass, double conc,
                              double omega_m) {
    if (radius <= 0 || mass <= 0 || conc <= 0) return std::numeric_limits<double>::quiet_NaN();
    const double rho_m0 = 2.775e11 * omega_m;
    const double r180 = std::cbrt(3.0 * mass / (4.0 * pi * 180.0 * rho_m0));
    const double rs = r180 / conc;
    const double factor = std::log1p(conc) - conc / (1.0 + conc);
    const double rho_s = mass / (4.0 * pi * rs*rs*rs * factor);
    const double norm = 2.0 * rs * rho_s;
    const double x = radius / rs;
    const double upper = std::min(x, 1.0);
    const double integral = gauss_legendre_integral(1e-6, upper) +
                            gauss_legendre_integral(1.0, x);
    const double sigma = norm * nfw_f(x);
    const double bar_sigma = 2.0 * norm * integral / (x*x);
    return (bar_sigma - sigma) / 1e12;
}

// ==========================================
// Function: assign_jackknife
// Method: Spherical k-means using BCG centers and deterministic seed.
// ==========================================
inline std::vector<int> assign_jackknife(const std::vector<Lens> &lenses, int n_jack) {
    const int n = static_cast<int>(lenses.size()), k = std::min(n_jack, n);
    std::vector<std::array<double, 3>> xyz(n), centers(k), next(k);
    for (int i = 0; i < n; ++i) xyz[i] = sky_xyz(lenses[i].ra, lenses[i].dec);
    std::mt19937 rng(42);
    std::vector<int> choices(n), labels(n);
    std::iota(choices.begin(), choices.end(), 0);
    std::shuffle(choices.begin(), choices.end(), rng);
    for (int j = 0; j < k; ++j) centers[j] = xyz[choices[j]];
    for (int iter = 0; iter < 100; ++iter) {
        std::vector<int> counts(k, 0);
        next.assign(k, {0, 0, 0});
        for (int i = 0; i < n; ++i) {
            double best = -2; int label = 0;
            for (int j = 0; j < k; ++j) {
                const double dot = xyz[i][0]*centers[j][0] + xyz[i][1]*centers[j][1] + xyz[i][2]*centers[j][2];
                if (dot > best) { best = dot; label = j; }
            }
            labels[i] = label; ++counts[label];
            for (int d = 0; d < 3; ++d) next[label][d] += xyz[i][d];
        }
        double max_change = 0;
        for (int j = 0; j < k; ++j) {
            if (counts[j] == 0) next[j] = xyz[choices[(j + iter + 1) % n]];
            const double norm = std::sqrt(next[j][0]*next[j][0] + next[j][1]*next[j][1] + next[j][2]*next[j][2]);
            for (int d = 0; d < 3; ++d) {
                next[j][d] /= norm;
                max_change = std::max(max_change, std::abs(next[j][d] - centers[j][d]));
            }
        }
        centers.swap(next);
        if (max_change < 1e-8) break;
    }
    for (int i = 0; i < n; ++i) {
        double best = -2;
        for (int j = 0; j < k; ++j) {
            const double dot = xyz[i][0]*centers[j][0] + xyz[i][1]*centers[j][1] + xyz[i][2]*centers[j][2];
            if (dot > best) { best = dot; labels[i] = j; }
        }
    }
    return labels;
}
} // namespace bcg
