"""Mixed coherent derivative of the SAME eleven-term conditional mean scalar.

``mixed`` is d²/(dλ dη), not half that derivative. Matter stages contain
physical per-cell Fourier amplitudes; Maxwell shift and momentum are GLOBAL.
The two formal directions are independent analytic directions (q,-q), so all
pairings are bilinear, NOT sesquilinear. No reality projection is performed.
This is neither ordered quantum covariance energy nor a finite-dose result.

For u_R = u0 + λ exp(iqR) u+ + η exp(-iqR) u- + λη ux,
C(q) acts on exp(iqR)u+. Consequently u- C(q) u+ and u+ C(-q) u-
are separate cross terms. C0 acts on mean/mixed. Onsite DeltaC is separate;
the matched physical C must appear once and the actual matched DeltaC is zero.
"""
import numpy as np


_STAGE_KEYS = frozenset(("mean", "plus", "minus", "mixed"))
_STAGE_FIELDS = ("density", "potential", "u", "p", "shift", "momentum",
                 "gradient", "field_gradient_electronic")


def _stages(stages, one_body, gamma):
    """Require actual complete inputs, including explicit non-Gamma zeros."""
    if set(stages) != _STAGE_KEYS:
        raise ValueError("stages require EXACT mean/plus/minus/mixed keys")
    if "mixed" not in one_body or one_body["mixed"] is None:
        raise ValueError("the actual mixed bare-electronic moment is required")
    if not isinstance(gamma, (bool, np.bool_)):
        raise TypeError("gamma must explicitly be boolean")
    for key in _STAGE_KEYS:
        for name in _STAGE_FIELDS:
            if not hasattr(stages[key], name) or getattr(stages[key], name) is None:
                raise ValueError(f"{key}.{name} must be supplied, including mixed data")
    if not gamma:
        for key in ("plus", "minus"):
            for name in ("shift", "momentum"):
                if np.any(np.asarray(getattr(stages[key], name)) != 0):
                    raise ValueError("non-Gamma first Maxwell directions must be explicit zero")
    return tuple(stages[key] for key in ("mean", "plus", "minus", "mixed"))


def _bilinear(left, right, pairing):
    """Product rule in order mean,plus,minus,mixed; no conjugation or half."""
    a, b, c, d = left
    e, f, g, h = right
    return pairing(a, h) + pairing(d, e) + pairing(b, g) + pairing(c, f)


def mixed_components(*, reference, coefficients, field, charges, stages,
                     one_body, grid_measure, curvature_plus, curvature_minus,
                     gamma):
    """Return complex-preserving mixed derivatives of all eleven components.

    Required stage attributes match ``joint_mean_energy_components.components``.
    ``one_body["mixed"]`` is the actual mixed bare-electronic scalar moment,
    never an electronic residual fitted to a total. Bare energy is linear in
    density: unused mean/first bare moments need not be computed. All four
    stages and both independent C(q)/C(-q) matrices are mandatory. ``gamma`` says
    whether first uniform-field directions are allowed; a non-Gamma mixed
    field remains allowed. Fixed coefficients, masses, charges, nbar and forces
    are NOT differentiated. No stationarity assumption removes mean×mixed.
    """
    r, c = reference, coefficients
    ss = _stages(stages, one_body, gamma)
    s0, sp, sm, sx = ss
    cp, cm = np.asarray(curvature_plus), np.asarray(curvature_minus)
    shape = np.asarray(r.curvature).shape
    if cp.shape != shape or cm.shape != shape:
        raise ValueError("independent C(q), C(-q) must have the C0 matrix shape")
    n = r.ncells
    get = lambda name: tuple(np.asarray(getattr(s, name)) for s in ss)
    dot = lambda a, b: a @ b
    mat = lambda matrix: lambda a, b: a @ matrix @ b
    density = list(get("density"))
    density[0] = density[0] - r.nbar
    za = np.repeat(charges, 3)
    kinetic = tuple(np.asarray(s.p) - za*np.tile(s.shift, len(charges)) for s in ss)
    u = get("u")
    return dict(
        bare_electronic=n*one_body["mixed"],
        hartree=n*grid_measure*_bilinear(density, get("potential"),
                                        lambda a, b: np.sum(a*b))/2,
        electron_ionic=n*_bilinear(u, get("gradient"), dot),
        ionic_reference_linear=-n*(sx.u @ (r.ionic_force+r.support_force)),
        ionic_C=n*(s0.u @ r.curvature @ sx.u + sx.u @ r.curvature @ s0.u
                   + sp.u @ cm @ sm.u + sm.u @ cp @ sp.u)/2,
        ionic_kinetic_minimal=n*_bilinear(kinetic, kinetic,
                                         lambda a, b: np.sum(a*b/r.masses))/2,
        delta_C=n*_bilinear(u, u, mat(c.delta_curvature))/2,
        electron_field=n*_bilinear(get("shift"), get("field_gradient_electronic"), dot),
        field_contact=n*_bilinear(get("shift"), get("shift"), mat(c.contact))/2,
        ion_field_mixed=n*_bilinear(u, get("shift"), mat(c.mixed)),
        maxwell=_bilinear(get("momentum"), get("momentum"), dot)/(2*field.mass))
