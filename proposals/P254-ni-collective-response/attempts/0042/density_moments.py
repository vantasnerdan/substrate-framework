"""All degree-14 finite-density moments without fifteen complex Bessel calls.

Same exact Legendre/Fourier identity as A31. Small arguments use the convergent
entire series; all remaining complex arguments descend from two independent
SciPy seeds. No unstable complex upward recurrence, gap clipping or real projection.
"""
import numpy as np
from scipy.special import spherical_jn

SMALL_SERIES_TERMS = 10


def spherical_orders(z, count):
    """Return j_0..j_(count-1) for literal complex arguments, bounded by order."""
    z = np.asarray(z, complex)
    reflect = z.real < 0
    w = np.where(reflect, -z, z).ravel()
    result = np.empty((count, w.size), complex)
    small = np.abs(w) <= 1
    regular = ~small
    if np.any(small):
        x = w[small]
        orders = np.arange(count)
        coefficients = np.ones((count, SMALL_SERIES_TERMS+1))
        for k in range(1, SMALL_SERIES_TERMS+1):
            coefficients[:, k] = coefficients[:, k-1]/(2*k*(2*orders+2*k+1))
        powers = np.empty((SMALL_SERIES_TERMS+1, len(x)), complex)
        powers[0] = 1
        square = -(x*x)
        for k in range(1, SMALL_SERIES_TERMS+1):
            powers[k] = powers[k-1]*square
        series = coefficients@powers
        leading = np.ones(len(x), complex)
        for order in range(count):
            result[order, small] = leading*series[order]
            leading *= x
            leading *= 1/(2*order+3)
    if np.any(regular):
        x = w[regular]
        current = spherical_jn(count-1, x)
        following = spherical_jn(count, x)
        result[count-1, regular] = current
        for order in range(count-1, 0, -1):
            previous = (2*order+1)*current/x-following
            result[order-1, regular] = previous
            following, current = current, previous
    for order in range(1, count, 2):
        result[order, reflect.ravel()] *= -1
    return result.reshape((count,)+z.shape)


def fourier_legendre_moments(gaps, half_width, midpoint, count):
    """Integral P_n((t-midpoint)/half_width)*exp(+i*gap*t), all n together."""
    gaps = np.asarray(gaps)
    moments = spherical_orders(gaps*half_width, count)
    phase_jacobian = half_width*np.exp(1j*gaps*midpoint)
    for order in range(count):
        moments[order] *= (2*(1j**order))*phase_jacobian
    return moments
