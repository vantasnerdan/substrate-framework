"""Conditional exact axisymmetric Euler/Grad--Shafranov material algebra.

The supplied streamfunction, toroidal label and pressure are treated as
actual symbolic fields on ``r>0``.  This module derives cylindrical velocity
and the literal Euler residual; it does not solve the finite-R inverse or
certify existence from placeholders.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class CompactRingFields:
    velocity: sp.ImmutableMatrix
    pressure: sp.Expr
    residual: sp.ImmutableMatrix


def compact_ring_fields(streamfunction, toroidal_label, pressure, r, z):
    """Return ``(u_r,u_phi,u_z)``, pressure and cylindrical Euler residual.

    Convention: ``u_r=-psi_z/r``, ``u_z=psi_r/r``, ``u_phi=I/r``.
    ``pressure`` is the supplied specific pressure (physical pressure divided
    by the constant density), in velocity-squared units. The residual is
    ``(u·grad)u + grad(pressure)`` including cylindrical metric terms. It
    vanishes for a supplied steady solution on ``r>0``. This function does
    not infer a Bernoulli label or derive the supplied pressure from one.
    Unknown symbolic signs/realness remain caller hypotheses.
    """
    r, z = sp.sympify(r), sp.sympify(z)
    if not isinstance(r, sp.Symbol) or not isinstance(z, sp.Symbol) or r == z:
        raise ValueError("r and z must be distinct symbolic coordinates")
    if r.is_nonpositive is True or r.is_real is False or z.is_real is False:
        raise ValueError("the cylindrical representation requires real coordinates and r > 0")
    psi, swirl_label, p = map(sp.sympify, (streamfunction, toroidal_label, pressure))
    ur, up, uz = -sp.diff(psi, z) / r, swirl_label / r, sp.diff(psi, r) / r
    rr = sp.simplify(ur*sp.diff(ur, r) + uz*sp.diff(ur, z) - up**2/r + sp.diff(p, r))
    rp = sp.simplify(ur*sp.diff(up, r) + uz*sp.diff(up, z) + ur*up/r)
    rz = sp.simplify(ur*sp.diff(uz, r) + uz*sp.diff(uz, z) + sp.diff(p, z))
    return CompactRingFields(sp.ImmutableMatrix([ur, up, uz]), p,
                              sp.ImmutableMatrix([rr, rp, rz]))


def isotropic_tag_normalization(velocity_square_integral, tag_second_moment,
                                tag_volume, cell_volume, density=1,
                                mean_tag_fraction=1):
    """Return literal ``(j,a,ell_tag_sq)`` geometric normalizations.

    ``velocity_square_integral`` is integral |u|^2 dx, without a rho/2
    factor: its units are length^5/time^2. ``tag_second_moment`` is
    integral chi |x-X|^2 dx (length^5), before the mean fraction is applied.
    Tag/cell volumes have length^3 units. Thus j has density*length^2
    units, a has velocity-squared units, and ell_tag_sq has length^2 units.

    Positive input integrals give positive normalizations; zero integrals
    give the corresponding degenerate value. Unknown symbolic signs remain
    caller hypotheses. These data do not assert a dynamical constitutive law.
    """
    speed_integral, moment, volume, cell, rho, frac = map(sp.sympify, (velocity_square_integral, tag_second_moment,
        tag_volume, cell_volume, density, mean_tag_fraction))
    for name, value in (("tag volume", volume), ("cell volume", cell)):
        if value.is_nonpositive is True or value.is_real is False:
            raise ValueError(f"{name} must be positive")
    if rho.is_nonpositive is True or rho.is_real is False:
        raise ValueError("density must be positive")
    if frac.is_nonpositive is True or (frac-1).is_positive is True or frac.is_real is False:
        raise ValueError("mean tag fraction must lie in (0,1]")
    if any(value.is_negative is True or value.is_real is False
           for value in (speed_integral, moment)):
        raise ValueError("squared-speed and tag-moment integrals must be nonnegative")
    return sp.ImmutableMatrix([
        2*rho*frac*moment/(3*cell), speed_integral/(3*cell), moment/volume
    ])
