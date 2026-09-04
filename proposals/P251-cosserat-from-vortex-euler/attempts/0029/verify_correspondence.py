"""Exact Euler operator and angle/rate correspondence audit, PR 199.

The Cartesian background is v0=Omega*(-y,x,0). Polar basis differentiation
derives both advective terms before solving for the mode velocity. No
historical residual operator is assumed. Units are tracked by independent
mass, length and time rescalings. Exact algebra only; no numerical threshold.
"""
import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0029-correspondence")
    r, rho, omega = s.symbols("r rho Omega", positive=True)
    theta, z, t = s.symbols("theta z t", real=True)
    m, k, w = s.symbols("m k w", real=True)
    wt = w - m * omega
    er = s.Matrix([s.cos(theta), s.sin(theta), 0])
    et = s.diff(er, theta)
    ez = s.Matrix([0, 0, 1])
    vr, vt, vz = s.symbols("v_r v_theta v_z")
    phase = s.exp(s.I * (m * theta + k * z - w * t))
    v = (vr * er + vt * et + vz * ez) * phase
    v0 = omega * r * et
    # (v0.grad)v = Omega*d_theta v; (v.grad)v0 includes basis derivatives.
    acceleration = (s.diff(v, t) + omega * s.diff(v, theta)
                    + phase * (vr * s.diff(v0, r)
                               + vt / r * s.diff(v0, theta)))
    components = s.Matrix([s.simplify(b.dot(acceleration) / phase)
                           for b in (er, et, ez)])
    expected = s.Matrix([-s.I * wt * vr - 2 * omega * vt,
                         -s.I * wt * vt + 2 * omega * vr,
                         -s.I * wt * vz])
    ledger.check("Cartesian Euler linearization fixes both Coriolis signs",
                 s.simplify(components - expected) == s.zeros(3, 1))
    p = s.Function("P")(r)
    gradient = s.Matrix([s.diff(p, r), s.I * m * p / r, s.I * k * p]) / rho
    solution = s.solve(components + gradient, [vr, vt, vz])
    D = 4 * omega**2 - wt**2
    original = s.Matrix([
        s.I * (wt * s.diff(p, r) - 2 * omega * m * p / r) / (rho * D),
        (2 * omega * s.diff(p, r) - wt * m * p / r) / (rho * D),
        k * p / (rho * wt),
    ])
    ledger.check("existing verify_cst002 fields solve the derived Euler operator",
                 s.simplify(original - s.Matrix([solution[q] for q in (vr, vt, vz)]))
                 == s.zeros(3, 1))
    historical = s.Matrix([
        s.I * (wt * s.diff(p, r) + 2 * omega * m * p / r) / (rho * D),
        -(2 * omega * s.diff(p, r) + wt * m * p / r) / (rho * D),
        k * p / (rho * wt),
    ])
    bad_residual = s.simplify((components + gradient).subs(
        dict(zip((vr, vt, vz), historical))))
    ledger.check("0019 reversed-Coriolis fields fail the actual Euler operator",
                 bad_residual != s.zeros(3, 1))
    print("0019 field residual:", bad_residual.T)

    # 0028 eta is the radial boundary displacement: eta has length units.
    Lv, a, eta, Om_i, Om_o, j = s.symbols("L_v a eta Omega_i Omega_o j", positive=True)
    M, L, T = s.symbols("mass_scale length_scale time_scale", positive=True)
    dimensions = {rho: M / L**3, Lv: L**-2, a: L, eta: L,
                  Om_i: T**-1, Om_o: T**-1, j: M / L}
    alpha_E = Lv * s.pi * rho * a**2 * eta**2 / 4
    alpha_gap = j * (Om_i - Om_o)**2 / 4
    energy_dim = s.simplify(alpha_E.subs(dimensions) / (s.pi / 4))
    # For the difference of rates, scale the entire expression, not each
    # rate to the same numerical value (which would manufacture a zero).
    rescaling = {q: q * dim for q, dim in dimensions.items()}
    gap_dim = s.simplify(alpha_gap.xreplace(rescaling) / alpha_gap)
    ledger.check("0028 alpha_E has microinertia units mass/length",
                 energy_dim == M / L)
    ledger.check("elastic alpha_gap has pressure units mass/(length*time^2)",
                 gap_dim == M / (L * T**2))
    ledger.check("equating the two alphas is dimensionally invalid",
                 s.simplify(energy_dim / gap_dim) == T**2)

    q = s.Function("relative_angle")(t)
    J, K = s.symbols("J K", positive=True)
    rate_energy = J * s.diff(q, t)**2 / 2
    angle_energy = K * q**2 / 2
    rate_el = s.diff(s.diff(rate_energy, s.diff(q, t)), t) - s.diff(rate_energy, q)
    angle_force = s.diff(angle_energy, q)
    ledger.check("rate-quadratic energy gives inertial acceleration",
                 s.simplify(rate_el - J * s.diff(q, t, 2)) == 0)
    ledger.check("angle-quadratic energy gives restoring torque",
                 s.simplify(angle_force - K * q) == 0)
    ledger.check("static rotation has no force from rate-quadratic energy",
                 rate_el.subs(q, s.Symbol("q0", real=True)).doit() == 0)

    # Rigid frame rotation in 3D leaves a Coriolis term. An axial inertial
    # plane wave has no pressure gradient transverse to its wavevector.
    wc = s.Symbol("relative_frequency", real=True)
    operator = s.Matrix([[-s.I * wc, -2 * Om_i], [2 * Om_i, -s.I * wc]])
    at_proposed_shift = s.factor(operator.det().subs(wc, 2 * (Om_i - Om_o)))
    ledger.check("finite-k background rotation is not a pure Doppler shift",
                 at_proposed_shift == 4 * Om_o * (2 * Om_i - Om_o))
    print("3D bulk residual of contrast-only inertial frequency:", at_proposed_shift)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
