"""Five-parameter group lensing fit called by the embedded Python runtime."""

import emcee
import numpy as np
from scipy.linalg import solve_triangular

from esd_nfw import OffsetNfw, nfw_esd
from esd_two_halo import TwoHalo, nonlinear_correlation
from esd_archive import PARAMETERS, save_fit
from esd_profiles import evaluate_profile_draws, stellar_profile

LOWER = np.array([12.0, 1.0, 0.0, 0.0, 0.3])
UPPER = np.array([16.0, 20.0, 1.0, 1.0, 10.0])


# ==========================================
# Class: HaloModel
# Method: Combine the fixed stellar mass with the mixed NFW and nonlinear two-halo terms.
# ==========================================
class HaloModel:
    # ==========================================
    # Function: __init__
    # Method: Size the universal NFW template for the allowed masses, radii and concentrations.
    # ==========================================
    def __init__(self, radius, stellar_mass, omega_m, correlation):
        self.stellar_mass = stellar_mass
        self.mean_density = 2.77536627e11 * omega_m
        min_x = np.min(radius) * LOWER[1]/self.halo_radius(10**UPPER[0])
        max_x = np.max(radius) * UPPER[1]/self.halo_radius(10**LOWER[0])
        self.offset = OffsetNfw(min_x, max_x, UPPER[1]*UPPER[3])
        self.correlation = correlation

    # ==========================================
    # Function: halo_radius
    # Method: Return comoving r180m in Mpc/h for a mass in M_sun/h.
    # ==========================================
    def halo_radius(self, mass):
        return np.cbrt(3*mass/(4*np.pi*180*self.mean_density))

    # ==========================================
    # Function: components
    # Method: Return one array per physical component for a batch of locations in parameter space.
    # ==========================================
    def components(self, theta, r, two_halo):
        log_mass, c, f_cen, eta, bias = np.atleast_2d(theta).T[:, :, None]
        mass = 10.0**log_mass
        r180 = self.halo_radius(mass)
        rs = r180/c
        normalization = mass / (2*np.pi*rs**2*(np.log1p(c) - c/(1+c))*1e12)
        x = r/rs
        centered = nfw_esd(x)
        off = self.offset.evaluate(x, eta*c, centered)
        stellar = stellar_profile(r, self.stellar_mass)
        one = normalization*(f_cen*centered + (1-f_cen)*off)
        two = bias * two_halo.evaluate(r180)
        return np.broadcast_to(stellar, one.shape), one, two


# ==========================================
# Class: FitProblem
# Method: Evaluate the joint Gaussian likelihood in whitened coordinates.
# ==========================================
class FitProblem:
    # ==========================================
    # Function: __init__
    # Method: Factor the full covariance once, and precompute the observed-radius model.
    # ==========================================
    def __init__(self, radius, measured, cov, stellar_mass, omega_m, corr_bytes):
        self.r = np.asarray(radius)
        self.y = np.asarray(measured)
        self.cov = np.asarray(cov).reshape((self.r.size, self.r.size))
        try:
            cholesky = np.linalg.cholesky(self.cov)
        except np.linalg.LinAlgError as exc:
            raise ValueError("The tangential jackknife covariance is not positive definite") from exc
        self.whitener = solve_triangular(cholesky, np.eye(self.r.size), lower=True)
        correlation = np.frombuffer(corr_bytes, dtype=np.float64).reshape(-1, 2)
        self.model = HaloModel(self.r, stellar_mass, omega_m, correlation)
        self.two_halo = TwoHalo(self.r, correlation, self.model.mean_density)

    # ==========================================
    # Function: log_probability
    # Method: Apply flat priors, then evaluate the complete non-diagonal likelihood.
    # ==========================================
    def log_probability(self, theta):
        theta = np.asarray(theta)
        in_prior = np.all((theta > LOWER) & (theta < UPPER), axis=1)
        log_p = np.full(theta.shape[0], -np.inf)
        if np.any(in_prior):
            model = sum(self.model.components(theta[in_prior], self.r, self.two_halo))
            residual = (self.y - model) @ self.whitener.T
            log_p[in_prior] = -0.5*np.sum(residual*residual, axis=1)
        return log_p


# ==========================================
# Function: sample_ensemble
# Method: Seed an independent emcee ensemble and return samples, log posterior, and diagnostics.
# ==========================================
def sample_ensemble(problem, walkers, steps, burn, seed, ensemble):
    rng = np.random.default_rng(np.random.SeedSequence([seed, ensemble]))
    sampler = emcee.EnsembleSampler(walkers, len(PARAMETERS), problem.log_probability, vectorize=True)
    sampler.random_state = np.random.RandomState(rng.integers(2**32)).get_state()
    initial = rng.uniform(LOWER, UPPER, size=(walkers, len(PARAMETERS)))
    sampler.run_mcmc(initial, steps, progress=False)
    chain = sampler.get_chain(discard=burn)
    logp = sampler.get_log_prob(discard=burn)
    if not np.all(np.isfinite(logp)):
        raise ValueError("Non-finite MCMC log posterior in retained chain")
    tau = (emcee.autocorr.integrated_time(chain, tol=0) if chain.shape[0] >= 8
           else np.full(len(PARAMETERS), np.nan))
    stats = np.r_[np.mean(sampler.acceptance_fraction), tau]
    packed = np.concatenate((chain, logp[:, :, None]), axis=2)
    return packed.astype(np.float64, copy=False).tobytes(), stats.astype(np.float64).tobytes()


# ==========================================
# Function: finish_fit
# Method: Capture the retained chain and joint prediction draws for a later pull.
# ==========================================
def finish_fit(problem, chain_bytes, stats_bytes, shape, archive, context):
    ensembles, steps, walkers = shape
    chain = np.frombuffer(chain_bytes, dtype=np.float64).reshape((ensembles, steps, walkers, 6))
    stats = np.frombuffer(stats_bytes, dtype=np.float64).reshape((ensembles, 6))
    flat = chain.reshape(-1, 6)
    best_index = int(np.argmax(flat[:, -1]))
    chi2 = float(-2.0 * flat[best_index, -1])
    dof = problem.r.size - len(PARAMETERS)
    curves = evaluate_profile_draws(problem, flat[:, :-1], best_index, context["profile_draws"], context["seed"])
    save_fit(archive, chain, stats, curves, context)
    print(f"Best retained draw: chi2/dof={chi2:.5g}/{dof}, samples={len(flat)}", flush=True)
    print(f"Saved fit archives: {archive}", flush=True)
