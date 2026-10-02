import numpy as np
from .validation import positive


def reflection(detuning, gamma_e, gamma_i):
    """Detuning in rad/s; gamma rates in s^-1, energy-linewidth convention."""
    positive('gamma_e', gamma_e, zero=True); positive('gamma_i', gamma_i, zero=True)
    d, e, i = np.broadcast_arrays(np.asarray(detuning), gamma_e, gamma_i)
    if not np.all(np.isfinite(d)):
        raise ValueError('Detuning must be finite.')
    numerator = 1j*d+(e-i)/2
    denominator = -1j*d+(e+i)/2
    # Fully decoupled lossless resonance: direct background r=-1, including d=0.
    out = np.full(d.shape, -1., dtype=complex)
    np.divide(numerator, denominator, out=out, where=denominator != 0)
    return out


def spectra(detuning, gamma_e, gamma_i):
    r = reflection(detuning, gamma_e, gamma_i)
    R = np.abs(r)**2
    return R, 1-R


def regime(gamma_e, gamma_i):
    positive('gamma_e', gamma_e, zero=True); positive('gamma_i', gamma_i, zero=True)
    if gamma_e == gamma_i == 0:
        return 'DECOUPLED / LOSSLESS (no absorption)'
    if np.isclose(gamma_e, gamma_i, rtol=1e-6, atol=1e-12):
        return 'CRITICAL COUPLING'
    return 'UNDER-COUPLED' if gamma_e < gamma_i else 'OVER-COUPLED'
