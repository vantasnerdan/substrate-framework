"""Canonical coordinate feedback for coherent physical mean variations.

First q/-q and mixed-Gamma electronic moments are supplied at their actual
variation order. Coordinates u,p are PHYSICAL per-cell amplitudes, not the
old collective sqrt(N) response variables. Maxwell a,Pi remain GLOBAL.
The fixed scalar is quadratic in canonical coordinates: all electronic
mixed forcing is already in the supplied mixed g and ell, not an extra
first-times-first force added here. This is not population noise or a
stationary approximation. C(q) and the original fixed DeltaC remain separate.
"""
import numpy as np


def coordinate_rhs(*, reference, coefficients, field, charges, polarizations,
                   curvature, gamma, small, gradient, field_gradient):
    """Return u,p,coordinate,momentum derivatives without a reality repair.

    The SAME map serves first and mixed variations; the mixed direction is
    Gamma. Non-Gamma first field coordinates/momenta must actually be empty.
    Source coefficient validation remains with the admitted material provider.
    No missing moment or field input is replaced by a guessed value.
    """
    masses = np.asarray(reference.masses)
    ndof = len(masses)
    u, p = np.asarray(small.u), np.asarray(small.p)
    a, pi = np.asarray(small.coordinate), np.asarray(small.momentum)
    g, ell = np.asarray(gradient), np.asarray(field_gradient)
    e = np.asarray(polarizations)
    if not isinstance(gamma, (bool, np.bool_)):
        raise TypeError('explicit Gamma selection required')
    nfield = len(e) if gamma else 0
    if (u.shape != (ndof,) or p.shape != (ndof,) or g.shape != (ndof,)
            or a.shape != (nfield,) or pi.shape != (nfield,) or ell.shape != (3,)):
        raise ValueError('complete physical variation coordinates/moments required')
    c = np.asarray(curvature) + coefficients.delta_curvature
    if c.shape != (ndof, ndof):
        raise ValueError('the actual C(q) matrix is required')
    force = -c @ u - g
    if not gamma:
        return dict(u=p/masses, p=force,
                    coordinate=np.zeros(0, complex), momentum=np.zeros(0, complex))
    charge = (np.repeat(charges, 3)[:, None]
              * np.tile(np.eye(3), (len(charges), 1))) @ e.T
    mixed = coefficients.mixed @ e.T
    contact = e @ coefficients.contact @ e.T
    velocity = (p-charge @ a)/masses
    force -= mixed @ a
    photon_force = -reference.ncells*(
        e @ ell + mixed.T @ u + contact @ a - charge.T @ velocity)
    return dict(u=velocity, p=force,
                coordinate=pi/field.mass, momentum=photon_force)
