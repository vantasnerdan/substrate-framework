"""Canonical total-force Fourier origins from DIRECT physical gradient data.

Input u and electronic gradient are physical per-cell; field coordinates are
global. sqrt(N) converts each force origin to the original collective Fourier
observable. No electronic origin is inferred as a residual of the total.
"""
import numpy as np


def force_components(curvature, mixed, u, coordinate, gradient, *, ncells):
    c, t = np.asarray(curvature), np.asarray(mixed)
    u, a, g = np.asarray(u), np.asarray(coordinate), np.asarray(gradient)
    if (u.ndim != 1 or g.shape != u.shape or c.shape != (len(u), len(u))
            or a.ndim != 1 or t.shape != (len(u), len(a))
            or not all(np.isfinite(x).all() for x in (c, t, u, a, g))
            or type(ncells) is not int or ncells < 1):
        raise ValueError('complete finite original force coordinates, matrices and positive N required')
    scale = np.sqrt(ncells)
    ionic, electronic, field = -c @ u, -g, -t @ a
    return dict(ionic=scale*ionic, electronic=scale*electronic,
                field=scale*field, total=scale*((ionic+electronic)+field))
