"""Physical BASELINE mean-energy ledger in global Ry, not heat or covariance.

This splits the existing JointReduction scalar without changing its model.
The minimal ionic kinetic replaces (does not supplement) neutral p^2/(2M).
Only Maxwell is already global; all matter/coupling terms carry ncells.
Array arithmetic deliberately preserves exact symbolic object scalars.
"""
import numpy as np


def components(*, reference, coefficients, field, charges, stage, one_body, grid_measure):
    """Return the eleven physical terms; no reality projection or fitted residual."""
    r, c, s = reference, coefficients, stage
    n = r.ncells
    za = np.repeat(charges, 3)*np.tile(s.shift, len(charges))
    return dict(
        bare_electronic=n*one_body,
        hartree=n*grid_measure*np.sum((s.density-r.nbar)*s.potential)/2,
        electron_ionic=n*(s.u @ s.gradient),
        ionic_reference_linear=-n*(s.u @ (r.ionic_force+r.support_force)),
        ionic_C=n*(s.u @ r.curvature @ s.u)/2,
        ionic_kinetic_minimal=n*np.sum((s.p-za)**2/r.masses)/2,
        delta_C=n*(s.u @ c.delta_curvature @ s.u)/2,
        electron_field=n*(s.shift @ s.field_gradient_electronic),
        field_contact=n*(s.shift @ c.contact @ s.shift)/2,
        ion_field_mixed=n*(s.u @ c.mixed @ s.shift),
        maxwell=(s.momentum @ s.momentum)/(2*field.mass))
