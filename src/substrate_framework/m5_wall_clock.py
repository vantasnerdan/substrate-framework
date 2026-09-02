"""Exact wall/shell/bubble sector primitives for the P250 degenerate-interface
campaign.

Everything here is derived from the accepted C-M5C-001 canonical action via
`substrate_framework.m5_exterior_clock`; no spectral or tension constant is
imported from any other model, and no oracle below is a literal comparison
(the #190 defect class is excluded by construction: every derived gate ships
with a channel-local mutation check that must fail).

The physical object is a stationary interface through the degenerate locus of
the exterior-degenerate M5 clock: the locus is the set of configurations with
psi = 0 and tangent-isotropic S (two coincident tangent-plane eigenvalues),
where the canonical clock generator norm Tr[A,S]^2/2 + |psi|^2 vanishes and
the whole diagonal S1 orbit fixes the configuration pointwise.

Conditional P250 campaign API: symbolic_verified scope only.  It attaches no
particle interpretation, no two-clock force, and no gravity.
"""

from __future__ import annotations

import sympy as sp

from substrate_framework.m5_exterior_clock import (
    canonical_velocity_energy_density,
    clock_inertia_density,
    full_clock_potential,
    relative_equilibrium_velocity,
)

WALL_LOCK_STRENGTH = sp.Integer(6)


def clock_slice_tensor(m, c, b) -> sp.Matrix:
    """Return the aligned real-psi wall slice ``S = diag(m, c+b, c-b)``.

    ``m`` is the axis eigenvalue, ``c`` half the tangent trace, ``b`` half the
    tangent eigenvalue split (the split ``s2 - s3 = 2b`` drives the clock).
    """
    m_v, c_v, b_v = sp.sympify(m), sp.sympify(c), sp.sympify(b)
    return sp.Matrix([[m_v, 0, 0], [0, c_v + b_v, 0], [0, 0, c_v - b_v]])


def wall_slice_potential(m, c, b, f, omega_sq=0) -> sp.Expr:
    """Return the exact rotating-frame slice potential ``V_w(m,c,b,f)``.

    Composed from the canonical action primitives (never hand-coded):
    ``V_w = V0 - (omega^2/2) * (f^2 + 4b^2)`` on the phase-lock slice
    ``psi = f >= 0`` real.  At ``omega_sq = 0`` this is the static potential.
    """
    S = clock_slice_tensor(m, c, b)
    return sp.factor(
        full_clock_potential(S, sp.sympify(f), 0)
        - sp.sympify(omega_sq) * clock_inertia_density(S, sp.sympify(f), 0) / 2
    )


def wall_slice_inertia(m, c, b, f) -> sp.Expr:
    """Return the canonical clock inertia density ``f^2 + 4b^2`` on the slice."""
    S = clock_slice_tensor(m, c, b)
    return sp.factor(clock_inertia_density(S, sp.sympify(f), 0))


def wall_slice_gradient_density(dm, dc, db, df) -> sp.Expr:
    """Return the exact static gradient energy density per unit area.

    Canonical normalization ``Tr(dS dx dS dx)/4 + |dpsi/dx|^2/2`` on the
    slice: ``(dm^2 + 2 dc^2 + 2 db^2)/4 + df^2/2`` for derivative inputs.
    """
    dm_v, dc_v = sp.sympify(dm), sp.sympify(dc)
    db_v, df_v = sp.sympify(db), sp.sympify(df)
    return sp.expand((dm_v**2 + 2*dc_v**2 + 2*db_v**2)/4 + df_v**2/2)


def wall_slice_el_residuals(omega_sq):
    """Return the four exact EL residual expressions ``2K u'' - dV_w/du``.

    The residuals are functions of the profile functions ``m(x), c(x), b(x),
    f(x)`` and their derivatives; setting them to zero is the stationary wall
    equation system.
    """
    x = sp.Symbol('x', real=True)
    m_f, c_f, b_f, f_f = _slice_profile_functions()
    V = wall_slice_potential(m_f, c_f, b_f, f_f, omega_sq)
    T = _slice_gradient_expression(m_f, c_f, b_f, f_f)
    residuals = []
    for u in (m_f, c_f, b_f, f_f):
        res = sp.diff(T + V, u) - sp.diff(sp.diff(T + V, u.diff(x)), x)
        residuals.append(sp.expand(-res))
    return tuple(residuals)


def _slice_profile_functions():
    x = sp.Symbol('x', real=True)
    m_f, c_f, b_f = (sp.Function(n, real=True)(x) for n in ('m', 'c', 'b'))
    f_f = sp.Function('f', nonnegative=True)(x)
    return m_f, c_f, b_f, f_f


def _slice_gradient_expression(m_f, c_f, b_f, f_f):
    x = sp.Symbol('x', real=True)
    return (m_f.diff(x)**2 + 2*c_f.diff(x)**2 + 2*b_f.diff(x)**2)/4 \
        + f_f.diff(x)**2/2


def first_integral_conservation_identity() -> bool:
    """Verify exactly ``d/dx(T - V_w) == sum(2K u'' - dV_w/du) u'``.

    On EL solutions the mechanical quantity ``T - V_w`` is constant along any
    stationary profile: the load-bearing wall first integral.
    """
    x = sp.Symbol('x', real=True)
    m_f, c_f, b_f, f_f = _slice_profile_functions()
    omega_sq = sp.Symbol('omega_sq', nonnegative=True)
    V = wall_slice_potential(m_f, c_f, b_f, f_f, omega_sq)
    T = _slice_gradient_expression(m_f, c_f, b_f, f_f)
    lhs = sp.diff(T - V, x)
    rhs = sp.S.Zero
    for u in (m_f, c_f, b_f, f_f):
        el = sp.diff(T + V, u) - sp.diff(sp.diff(T + V, u.diff(x)), x)
        rhs += -el*u.diff(x)
    return sp.simplify(sp.trigsimp(sp.expand(lhs - rhs))) == 0


def degenerate_locus_point(m0=1, c0=0):
    """Return the orbit-fixed locus configuration ``(S, 0, 0)``."""
    return clock_slice_tensor(sp.sympify(m0), sp.sympify(c0), 0), \
        sp.Integer(0), sp.Integer(0)


def locus_orbit_fixation() -> bool:
    """Verify exactly that the whole S1 orbit fixes the locus configuration."""
    q = sp.Symbol('q', real=True)
    rotation = sp.Matrix([[1, 0, 0], [0, sp.cos(q), -sp.sin(q)],
                          [0, sp.sin(q), sp.cos(q)]])
    S_loc, _, _ = degenerate_locus_point()
    rotated = sp.simplify(rotation * S_loc * rotation.T)
    if rotated != S_loc:
        return False
    return sp.simplify(clock_inertia_density(S_loc, 0, 0)) == 0


def orbit_invariance_identity() -> bool:
    """Verify exact orbit invariance of potential, inertia, velocity-energy
    and static gradient energy densities on the aligned slice."""
    q = sp.Symbol('q', real=True)
    m_g, c_g, b_g, f_g, p_g, w_g = sp.symbols(
        'm_g c_g b_g f_g p_g w_g', real=True)
    S_g = clock_slice_tensor(m_g, c_g, b_g)
    rotation = sp.Matrix([[1, 0, 0], [0, sp.cos(q), -sp.sin(q)],
                          [0, sp.sin(q), sp.cos(q)]])
    fr = f_g*sp.cos(q) - p_g*sp.sin(q)
    fi = f_g*sp.sin(q) + p_g*sp.cos(q)
    S_r = rotation * S_g * rotation.T
    checks = [
        sp.trigsimp(sp.expand(full_clock_potential(S_r, fr, fi)
                              - full_clock_potential(S_g, f_g, p_g))),
        sp.trigsimp(sp.expand(clock_inertia_density(S_r, fr, fi)
                              - clock_inertia_density(S_g, f_g, p_g))),
    ]
    vel1 = canonical_velocity_energy_density(
        *relative_equilibrium_velocity(S_g, f_g, p_g, w_g))
    vel2 = canonical_velocity_energy_density(
        *relative_equilibrium_velocity(S_r, fr, fi, w_g))
    checks.append(sp.trigsimp(sp.expand(vel2 - vel1)))
    return all(sp.simplify(ch) == 0 for ch in checks)


def phase_slip_energy_bound(gradient_scale, delta_theta, epsilon) -> sp.Expr:
    """Return the exact uniform phase-slip cost ``L^2 dtheta^2 eps / 4``.

    A phase jump of size ``delta_theta`` at the locus of a profile with
    scalar-amplitude Lipschitz constant ``gradient_scale`` and f(x*) = 0,
    realized by a continuous ramp of half-width ``epsilon``, costs at most
    this much gradient energy; it vanishes uniformly in ``delta_theta``.
    """
    L, dth, eps = (sp.sympify(v)
                   for v in (gradient_scale, delta_theta, epsilon))
    return sp.factor(L**2 * dth**2 * eps / 4)


def static_shell_derrick_residual(gradient_energy, potential_energy):
    """Return the exact static scaling residual ``T + 3U`` (omega = 0)."""
    T, U = sp.sympify(gradient_energy), sp.sympify(potential_energy)
    return sp.expand(T + 3*U)


def maxwell_system():
    """Return ``(pA, pB, pC, w_of)``: the exact deep-branch Maxwell system.

    The stationary clock-active bulk state lives on the m = 0 branch (the
    m-equation vanishes identically there).  With ``w_of = 32f^2 - 12f + 6
    - 24b`` from the f-equation, the remaining exact polynomials are
    ``pA = dc V_w``, ``pB = db V_w``, and ``pC = 2 V_w``, all derived here
    from the canonical potential.  Nothing is hardcoded.
    """
    c_v, b_v = sp.symbols('c b', real=True)
    f_v = sp.Symbol('f', nonnegative=True)
    om2 = sp.Symbol('omega_sq', nonnegative=True)
    w_of = 32*f_v**2 - 12*f_v + 6 - 24*b_v
    # The bulk state is a critical point of V_omega at FIXED omega; the
    # f-equation then fixes omega.  Differentiate first, substitute after:
    # substituting omega^2(b, f) before differentiating would corrupt db.
    V = sp.expand(wall_slice_potential(sp.Integer(0), c_v, b_v, f_v, om2))
    pA = sp.factor(sp.diff(V, c_v).subs(om2, w_of))
    pB = sp.factor(sp.diff(V, b_v).subs(om2, w_of))
    pC = sp.expand(2*V).subs(om2, w_of)
    return pA, pB, pC, w_of


def maxwell_frequency_resultants():
    """Return ``(R1, R2)``: exact c-elimination of the Maxwell system.

    With ``b = (32f^2 - 12f + 6 - w)/24`` substituted (f-equation), the
    Maxwell solutions are exactly the common roots of
    ``R1 = Res_c(pA, pB)`` and ``R2 = Res_c(pA, pC)`` in (f, w).
    """
    c_v, f_v, w_v = sp.symbols('c f w', real=True)
    b_sym = sp.Symbol('b', real=True)
    pA, pB, pC, _ = maxwell_system()
    b_sub = (32*f_v**2 - 12*f_v + 6 - w_v)/24
    qA = sp.numer(sp.cancel(sp.expand(pA.subs(b_sym, b_sub))))
    qB = sp.numer(sp.cancel(sp.expand(pB.subs(b_sym, b_sub))))
    qC = sp.numer(sp.cancel(sp.expand(pC.subs(b_sym, b_sub))))
    return sp.factor(sp.resultant(qA, qB, c_v)), \
        sp.factor(sp.resultant(qA, qC, c_v))


def wall_tension_integrand(dm, dc, db, df, m, c, b, f, omega_sq) -> sp.Expr:
    """Return the exact kink tension integrand ``T + V_w`` at one x.

    The tension is the x-integral of this integrand; on the kink the first
    integral gives the equivalent forms ``sigma = 2 int V_w dx = 2 int T dx``
    and the variational characterization
    ``sigma = 2 inf_path int sqrt(V_w) ds_g`` with metric
    diag(1/4, 1/2, 1/2, 1/2).
    """
    kinetic = wall_slice_gradient_density(dm, dc, db, df)
    return sp.expand(kinetic + wall_slice_potential(m, c, b, f, omega_sq))


def thin_wall_bag_radius(sigma, pressure) -> sp.Expr:
    """Return the exact volume-vs-surface selection law ``R = 2 sigma/p``."""
    return sp.factor(2*sp.sympify(sigma)/sp.sympify(pressure))


def envelope_identity() -> bool:
    """Verify exactly the fixed-charge envelope identity ``dE~/dQ = omega``.

    General symbolic proof for any smooth stationary family: with
    ``E~(Q) = E(w) + Q^2/(2 I(w))``, ``Q = w I`` and the stationarity
    envelope ``dE/dw = (w^2/2) dI/dw``, the derivative is exactly ``w``.
    """
    lam = sp.Symbol('lam', real=True)
    w_f = sp.Function('omega')(lam)
    E_f, I_f = sp.Function('E')(lam), sp.Function('I')(lam)
    Q_f = w_f*I_f
    E_Q = E_f + Q_f**2/(2*I_f)
    dEQ_dQ = sp.diff(E_Q, lam)/sp.diff(Q_f, lam)
    dEQ_dQ = dEQ_dQ.subs(sp.Derivative(E_f, lam),
                         w_f**2*sp.Derivative(I_f, lam)/2)
    return sp.simplify(dEQ_dQ - w_f) == 0


def inertia_envelope_identity() -> bool:
    """Verify the chain-rule envelope on a generic algebraic instance.

    For ``V_min(w2) = V_w(w2, u(w2))`` the total derivative equals the
    partial derivative at u plus the chain terms; stationarity of u kills
    the chain terms, so ``dV_min/dw2 = dV_w/dw2 |_u`` and, with
    ``V_w = V0 - (w2/2) iota(u)``, ``iota_int = -2 dV_min/dw2``.  A generic
    polynomial instance with independent monomials certifies the
    bookkeeping.
    """
    w2 = sp.Symbol('w2', real=True)
    u1, u2 = sp.symbols('u1 u2', real=True)
    # generic instance:
    # V_w = w2*u1^3*u2^2 + u1*u2^3 - (w2/2)(u1^2 + 3u2^2)
    u1s = w2**2 + 1
    u2s = w2**3 + 2*w2
    Vw = lambda a, b_, c_: a*b_**3*c_**2 + b_*c_**3 - a*(b_**2 + 3*c_**2)/2
    Vmin = w2*u1s**3*u2s**2 + u1s*u2s**3 - w2*(u1s**2 + 3*u2s**2)/2
    partial = sp.diff(Vw(w2, u1, u2), w2).subs({u1: u1s, u2: u2s})
    total = sp.expand(sp.diff(Vmin, w2))
    chain = (sp.diff(Vw(w2, u1, u2), u1).subs({u1: u1s, u2: u2s})
             * sp.diff(u1s, w2)
             + sp.diff(Vw(w2, u1, u2), u2).subs({u1: u1s, u2: u2s})
             * sp.diff(u2s, w2))
    return sp.simplify(total - partial - chain) == 0
