"""Nonlinear matter correlation and the projected two-halo term."""

import numpy as np
from scipy.fft import fht, fhtoffset
from scipy.integrate import cumulative_simpson
from scipy.interpolate import PchipInterpolator


# ==========================================
# Function: nonlinear_correlation
# Method: Normalize CAMB to sigma8 and transform the original-Halofit spectrum with FFTLog.
# ==========================================
def nonlinear_correlation(omega_m, h, omega_b_h2, spectral_index, sigma8, neutrino_mass, z):
    import camb
    from camb import model

    p = camb.CAMBparams()
    cold = omega_m*h*h - omega_b_h2 - neutrino_mass/93.14
    p.set_cosmology(H0=100*h, ombh2=omega_b_h2, omch2=cold, mnu=neutrino_mass, omk=0)
    p.InitPower.set_params(As=2.1e-9, ns=spectral_index)
    p.WantCls = False
    p.set_matter_power(redshifts=[z, 0.0], kmax=100*h, nonlinear=True)
    p.NonLinear = model.NonLinear_pk
    p.NonLinearModel.set_params(halofit_version="original")
    result = camb.get_results(p)
    amplitude = 2.1e-9 * (sigma8 / result.get_sigma8_0())**2
    p.InitPower.set_params(As=amplitude, ns=spectral_index)
    result.power_spectra_from_transfer(p.InitPower)
    power = result.get_matter_power_interpolator(
        nonlinear=True, hubble_units=True, k_hunit=True, extrap_kmax=1000
    )
    k = np.geomspace(max(power.kmin, 1e-5), 1000, 4096)
    dlog = np.log(k[1]/k[0])
    offset = fhtoffset(dlog, mu=0.5)
    r = np.exp(offset) / k[::-1]
    # Suppress endpoint wrap-around outside the scales used by the projection.
    log_k = np.log(k)
    window_lo = np.clip((log_k - log_k[0]) / np.log(8), 0, 1)
    window_hi = np.clip((log_k[-1] - log_k) / np.log(10), 0, 1)
    window = np.sin(np.pi*window_lo/2)**2 * np.sin(np.pi*window_hi/2)**2
    transformed = fht(k**1.5 * power.P(z, k) * window, dlog, mu=0.5, offset=offset)
    xi = np.sqrt(np.pi/2) * transformed / (2*np.pi**2*r**1.5)
    use = (r >= 0.018) & (r <= 1200)
    r, xi = r[use], xi[use]
    zeta = (1 + 1.17*xi)**1.49 / (1 + 0.69*xi)**2.09
    threshold = 10.0**0.9  # Equation A8, using a base-10 log.
    crossings = np.flatnonzero((zeta[:-1]*xi[:-1] >= threshold) & (zeta[1:]*xi[1:] < threshold))
    if crossings.size == 0 or not np.all(np.isfinite(zeta)):
        raise ValueError("CAMB correlation does not bracket the scale-dependent bias transition")
    i = crossings[-1]
    frozen = np.interp(threshold, (zeta[i+1]*xi[i+1], zeta[i]*xi[i]), (zeta[i+1], zeta[i]))
    zeta[:i+1] = frozen
    return np.ascontiguousarray(np.column_stack((r, zeta*xi)), dtype=np.float64).tobytes()


# ==========================================
# Class: TwoHalo
# Method: Project the biased matter correlation, excluding a sphere around the host.
# ==========================================
class TwoHalo:
    # ==========================================
    # Function: __init__
    # Method: Tabulate the projection in u=asinh(LOS/R) and integrate the inner spherical term.
    # ==========================================
    def __init__(self, radius, correlation, mean_density):
        r, biased = correlation.T
        if np.max(radius) >= r[-1]:
            raise ValueError("Fit radii exceed the two-halo projection grid")
        self.radius = radius
        self.u_max = np.arccosh(r[-1]/radius)
        self.scale = mean_density / 1e12
        self.points = 513
        u = self.u_max[:, None] * np.linspace(0, 1, self.points)
        sphere = PchipInterpolator(r, biased)
        rho = sphere(np.clip(radius[:, None]*np.cosh(u), r[0], r[-1]))
        integrand = rho * radius[:, None] * np.cosh(u) * np.exp(-2*u)
        self.tail = cumulative_simpson(integrand[:, ::-1], x=-u[:, ::-1],
                                       axis=1, initial=0)[:, ::-1]
        self.primitive = PchipInterpolator(r, cumulative_simpson(biased*r*r, x=r, initial=0))
        self.inside = self.primitive(np.clip(radius, r[0], r[-1]))

    # ==========================================
    # Function: evaluate
    # Method: Combine the projected shell and spherical integrals at a mass-dependent cutoff.
    # ==========================================
    def evaluate(self, cutoff):
        r = self.radius
        u_min = np.arccosh(np.maximum(1, cutoff/r))
        location = np.clip(u_min/self.u_max*(self.points-1), 0, self.points-1)
        lo = np.minimum(location.astype(int), self.points-2)
        fraction = location - lo
        radial = np.arange(r.size)
        tail = (1-fraction)*self.tail[radial, lo] + fraction*self.tail[radial, lo+1]
        inside = np.where(r > cutoff, 4*(self.inside - self.primitive(cutoff))/(r*r), 0)
        return self.scale * (inside - 2*tail)
