import numpy as np
from scipy.constants import epsilon_0, mu_0
from .validation import positive, integer


def cutoff(w, h, eps_r=1., m=1, n=0):
    for name, value in [('w', w), ('h', h), ('eps_r', eps_r)]:
        positive(name, value)
    integer('m', m); integer('n', n)
    if m == n == 0:
        raise ValueError('Rectangular PEC guide has no TE00 channel.')
    return np.hypot(m/w, n/h) / (2*np.sqrt(mu_0*epsilon_0*eps_r))


def beta10(f, w, h, eps_r=1.):
    positive('frequency', f)
    fc = cutoff(w, h, eps_r)
    ratio = np.asarray(f)/fc
    # Stable at cutoff; positive imaginary branch for evanescence.
    return np.pi/w * np.lib.scimath.sqrt((ratio-1)*(ratio+1))


def guide_wavelength(beta):
    if np.any(np.imag(beta) != 0) or np.any(np.real(beta) <= 0):
        raise ValueError('Guide wavelength requires a propagating channel above cutoff.')
    return 2*np.pi/np.real(beta)


def single_te10(f, w, h, eps_r=1.):
    # Every other TE/TM cutoff is >= min(TE20, TE01).
    return cutoff(w,h,eps_r) < f < min(cutoff(w,h,eps_r,2,0), cutoff(w,h,eps_r,0,1))
