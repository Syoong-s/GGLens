"""Centered and Rayleigh-miscentered NFW profiles in dimensionless radii."""

import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import RectBivariateSpline
from scipy.special import i0e, roots_legendre


# ==========================================
# Function: nfw_projection
# Method: Evaluate the analytic projection and its regular limit at x=1.
# ==========================================
def nfw_projection(x):
    x = np.asarray(x, dtype=float)
    t = x*x - 1.0
    near = np.abs(t) < 1e-4
    inner = (t < 0) & ~near
    outer = (t > 0) & ~near
    h = np.empty_like(x)
    h[near] = 1 + t[near]*(-1/3 + t[near]*(1/5 - t[near]/7))
    h[inner] = np.arccosh(1.0/x[inner]) / np.sqrt(-t[inner])
    h[outer] = np.arctan(np.sqrt(t[outer])) / np.sqrt(t[outer])
    surface = np.empty_like(x)
    surface[near] = 1/3 + t[near]*(-1/5 + t[near]/7)
    surface[~near] = (1.0 - h[~near]) / t[~near]
    return surface, h


# ==========================================
# Function: nfw_esd
# Method: Subtract the analytic local surface density from its interior mean.
# ==========================================
def nfw_esd(x):
    x = np.asarray(x, dtype=float)
    surface, h = nfw_projection(x)
    small = x < 0.08
    g = np.array(np.log(x / 2.0) + h)
    z = x[small]
    s = np.sqrt(1.0 - z*z)
    # Rearrange the two large logarithms before dividing by x^2.
    g[small] = z*z * np.log(2.0/z) / (s*(1+s)) + np.log1p(-z*z / (2*(1+s))) / s
    return 2.0 * g / (x*x) - surface


# ==========================================
# Function: offset_surface
# Method: Integrate the NFW surface density convolved with a 2D isotropic Gaussian.
# ==========================================
def offset_surface(x, sigma, nodes, weights):
    lo = np.maximum(-x/sigma, -9.0)
    half = (9.0 - lo) / 2.0
    t = lo[:, None] + half[:, None] * (nodes + 1)
    u = x[:, None] + sigma * t
    # The exponential in I0 cancels the linear cross term in the Gaussian.
    kernel = np.exp(-0.5*t*t) * i0e(u*x[:, None]/(sigma*sigma)) * u/sigma
    return half * np.sum(weights * kernel * nfw_projection(u)[0], axis=1)


# ==========================================
# Class: OffsetNfw
# Method: Tabulate excess density after Rayleigh averaging of the projected density.
# ==========================================
class OffsetNfw:
    # ==========================================
    # Function: __init__
    # Method: Integrate a log grid after subtracting the constant central surface density.
    # ==========================================
    def __init__(self, min_x, max_x, max_sigma):
        self.min_sigma = min(0.005, min_x/20)
        log_sigma = np.linspace(np.log(self.min_sigma), np.log(max_sigma), 128)
        log_x = np.linspace(np.log(min_x/256), np.log(max_x*1.01), 640)
        x = np.exp(log_x)
        nodes, weights = roots_legendre(96)
        table = np.empty((log_sigma.size, x.size))
        for i, sigma in enumerate(np.exp(log_sigma)):
            baseline = offset_surface(np.array([0.0]), sigma, nodes, weights)[0]
            surface = offset_surface(x, sigma, nodes, weights)
            difference = surface - baseline
            inner = x[0]**2 * difference[0] / 4.0
            diff_area = inner + cumulative_simpson(difference*x*x, x=log_x, initial=0)
            area = inner + baseline*x[0]**2/2.0
            area += cumulative_simpson(surface*x*x, x=log_x, initial=0)
            # Removing the constant is stable near the center; the direct integral
            # avoids subtracting two large areas in the outer halo.
            table[i] = np.where(x < sigma/4, 2*diff_area/(x*x) - difference,
                                2*area/(x*x) - surface)
        self.table = RectBivariateSpline(log_sigma, log_x, table)

    # ==========================================
    # Function: evaluate
    # Method: Interpolate the offset profile, taking the centered limit for vanishing width.
    # ==========================================
    def evaluate(self, x, sigma, centered):
        sigma = np.broadcast_to(sigma, x.shape)
        result = centered.copy()
        valid = sigma >= self.min_sigma
        result[valid] = self.table.ev(np.log(sigma[valid]), np.log(x[valid]))
        return result
