"""Explicit display-to-SI conversions."""
import numpy as np

def ghz_to_hz(x): return np.asarray(x) * 1e9

def mm_to_m(x): return np.asarray(x) * 1e-3

def hz_to_ghz(x): return np.asarray(x) / 1e9

def m_to_mm(x): return np.asarray(x) * 1e3

def hz_to_omega(x): return 2 * np.pi * np.asarray(x)
