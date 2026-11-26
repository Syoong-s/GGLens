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

// East and north share a cos(source Dec) factor, which cancels in the bearing.
inline double source_bearing(const double *source, const std::array<double, 3> &lens,
                             double xy, double east) {
    const double north = (source[0] * source[0] + source[1] * source[1]) * lens[2]
        - source[2] * xy;
    return std::atan2(north, east);
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
    static double beta(double zl, double zs, double chil, double chis) {
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

    // Radius queries emit a subsequence of this order, including after cuts.
    int preorder_point(std::size_t i) const { return nodes_[i].point; }

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

} // namespace bcg
