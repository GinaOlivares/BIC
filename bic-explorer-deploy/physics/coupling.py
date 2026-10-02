import numpy as np
from .validation import positive, integer


def normalized_coupling(beta, length, m):
    integer('m', m); positive('L', length, zero=True)
    if np.any(np.imag(beta) != 0):
        raise ValueError('Wall-interference model requires real propagating beta.')
    positive('beta', beta)
    return 1 - (-1)**m*np.exp(2j*np.asarray(beta)*np.asarray(length))


def radiation(beta, length, m, gamma_plus=1.):
    positive('gamma_plus', gamma_plus)
    return gamma_plus * np.abs(normalized_coupling(beta, length, m))**2
