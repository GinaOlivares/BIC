import numpy as np


def positive(name, value, zero=False):
    x = np.asarray(value)
    if not np.all(np.isfinite(x)) or np.any(x < 0 if zero else x <= 0):
        raise ValueError(f'{name} must be finite and {"nonnegative" if zero else "positive"}.')


def integer(name, value):
    positive(name, value, zero=True)
    if int(value) != value:
        raise ValueError(f'{name} must be an integer.')


def geometry(w, h, a, b):
    for name, value in [('w', w), ('h', h), ('a', a), ('b', b)]:
        positive(name, value)
    if a >= b or 2*b >= w:
        raise ValueError('Require 0 < a < b < w/2.')
