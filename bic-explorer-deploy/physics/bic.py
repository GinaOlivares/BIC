import numpy as np
from .validation import positive, integer
from .waveguide import guide_wavelength

BIC_TOL = 1e-12  # normalized squared coupling; a plotting/diagnostic threshold


def position(beta, m, q):
    integer('m', m); integer('q', q)
    return guide_wavelength(beta)*(q/2 + (m % 2)/4)


def quadratic_rate(beta, length, lb, gamma_plus=1.):
    positive('gamma_plus', gamma_plus); positive('beta', beta)
    return 4*gamma_plus*beta**2*(np.asarray(length)-lb)**2


def pa_positions(beta, lb, gamma_i, gamma_plus=1., exact=False):
    positive('gamma_i', gamma_i, zero=True); positive('gamma_plus', gamma_plus)
    positive('beta', beta)
    if gamma_i == 0 or gamma_i > 4*gamma_plus:
        return []  # zero loss cannot absorb; no critical point above maximum rate
    v = .5*np.sqrt(gamma_i/gamma_plus)
    delta = (np.arcsin(v) if exact else v)/beta
    return sorted(set(float(x) for x in [lb-delta, lb+delta] if x >= 0))


def log_quality(omega0, gamma_e, gamma_plus=1., cap=16.):
    positive('omega0', omega0); positive('gamma_e', gamma_e, zero=True)
    positive('gamma_plus', gamma_plus)
    rate = np.asarray(gamma_e, dtype=float)
    mask = rate/gamma_plus <= BIC_TOL
    # Log domain avoids overflow. Masked points are explicitly marked by the UI.
    logq = np.log10(omega0)-np.log10(np.maximum(rate, np.finfo(float).tiny))
    return np.where(mask, cap, np.minimum(logq, cap)), mask
