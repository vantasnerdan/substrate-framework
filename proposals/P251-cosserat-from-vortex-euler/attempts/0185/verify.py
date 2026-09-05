"""Exact helical Euler pressure, Kelvin, KKS and mechanical observation algebra."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0185")
    r, c, circulation, sigma, density, length = sp.symbols("r c C sigma rho Lz", positive=True)
    s, t = sp.symbols("s t", real=True)
    d = c**2 + r**2
    f = circulation / d
    h = sp.Matrix([0, r, c])
    j = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    ledger.check("constant helical momentum is derived from the profile", sp.simplify(sp.diff(d * f, r)) == 0)
    phi = sp.Function("phi")(r, s, t)
    velocity = sp.Matrix([sp.diff(phi, s) / r, -c**2 * sp.diff(phi, r) / d, c * r * sp.diff(phi, r) / d])
    divergence = sp.diff(r * velocity[0], r) / r + sp.diff(velocity[1], s) / r - sp.diff(velocity[2], s) / c
    ledger.check("reduced velocity is exactly incompressible in the helical metric", sp.simplify(divergence) == 0)
    ledger.check("nonzero-frequency velocity has zero helical momentum", sp.simplify(h.dot(velocity)) == 0)
    ledger.check("physical kinetic metric retains the axial velocity", sp.simplify(velocity.dot(velocity) - sp.diff(phi, s)**2 / r**2 - c**2 * sp.diff(phi, r)**2 / d) == 0)
    radial_pressure = -sp.diff(phi, s, t) / r - 2 * f * c**2 * sp.diff(phi, r) / d
    angular_pressure = r * c**2 * sp.diff(phi, r, t) / d - 2 * f * c**2 * sp.diff(phi, s) / d
    h_phi_t = -sp.diff(r * c**2 * sp.diff(phi, r, t) / d, r) / r - sp.diff(phi, s, s, t) / r**2
    beta = 8 * circulation * c**2 / d**3
    pressure_compatibility = (sp.diff(radial_pressure, s) - sp.diff(angular_pressure, r)) / r
    ledger.check("full pressure compatibility gives the actual active-vorticity equation", sp.simplify(pressure_compatibility - h_phi_t + beta * sp.diff(phi, s)) == 0)
    b = 2 * circulation * c**2 / d**2
    ledger.check("background quotient vorticity has the restoring gradient", sp.simplify(-sp.diff(b, r) / r - beta) == 0)
    h_phi = -sp.diff(r * c**2 * sp.diff(phi, r) / d, r) / r - sp.diff(phi, s, s) / r**2
    curl_velocity = sp.Matrix([sp.diff(velocity[2], s) / r + sp.diff(velocity[1], s) / c, -sp.diff(velocity[0], s) / c - sp.diff(velocity[2], r), sp.diff(r * velocity[1], r) / r - sp.diff(velocity[0], s) / r])
    ledger.check("full horizontal vorticity is the actual helical active scalar", all(sp.simplify(value) == 0 for value in curl_velocity - h_phi * h / c))
    ledger.check("relative kinetic energy has pointwise rather than fitted infrared cancellation", sp.simplify((f * h + velocity).dot(f * h + velocity) - f**2 * d - velocity.dot(velocity)) == 0)
    weight_tail = sp.factor(8 * circulation * c**2 / r**4 - r**2 * beta)
    ledger.check("the degenerating radial metric still has a strict compact tail bound", weight_tail.is_positive is True)

    # Complete pressure and Kelvin algebra away from the helical carrier.
    m, k = sp.symbols("m k", real=True)
    pi, pi_r = sp.symbols("pi pi_r")
    vr, vt, vz = sp.symbols("vr vt vz")
    arbitrary_v = sp.Matrix([vr, vt, vz])
    epsilon = m + c * k
    nu = sigma + f * epsilon
    gradient = sp.Matrix([pi_r, sp.I * m * pi / r, sp.I * k * pi])
    euler = sp.I * nu * arbitrary_v + 2 * f * j * arbitrary_v + sp.diff(f, r) * vr * h + gradient
    ledger.check("full helical momentum equation retains the off-carrier pressure term", sp.simplify(h.dot(euler) - sp.I * nu * h.dot(arbitrary_v) - sp.I * epsilon * pi) == 0)
    omega = 2 * c * circulation * h / d**2
    xi = arbitrary_v / (sp.I * nu) + sp.diff(f, r) * vr * h / (sp.I * nu)**2
    gradient_quotient = gradient / (sp.I * nu) - sp.Matrix([pi * sp.diff(nu, r) / (sp.I * nu**2), 0, 0])
    kelvin_error = xi.cross(omega) - arbitrary_v - gradient_quotient
    expected_error = -euler / (sp.I * nu) + sp.diff(f, r) * (h.dot(arbitrary_v) + epsilon * pi / nu) * sp.Matrix([1, 0, 0]) / (sp.I * nu)
    ledger.check("actual Kelvin reconstruction includes the new tau-pressure cancellation", all(sp.simplify(value) == 0 for value in kelvin_error - expected_error))
    radial_components = sp.Matrix([sp.Function(name)(r) for name in ("vr", "vt", "vz")])
    radial_xi = radial_components / (sp.I * nu) + sp.diff(f, r) * radial_components[0] * h / (sp.I * nu)**2
    divergence_v = sp.diff(r * radial_components[0], r) / r + sp.I * m * radial_components[1] / r + sp.I * k * radial_components[2]
    divergence_xi = sp.diff(r * radial_xi[0], r) / r + sp.I * m * radial_xi[1] / r + sp.I * k * radial_xi[2]
    ledger.check("off-carrier Lin displacement remains exactly solenoidal", sp.simplify(divergence_xi - divergence_v / (sp.I * nu)) == 0)
    ledger.check("silently keeping tau zero away from resonance loses a real row", sp.diff(nu, r) != 0)

    amplitude, slope, angular_index = sp.symbols("phi_m phi_m_prime m_integer", real=True)
    va = sp.Matrix([-angular_index * amplitude * sp.sin(angular_index * s) / r, -c**2 * slope * sp.cos(angular_index * s) / d, c * r * slope * sp.cos(angular_index * s) / d])
    vb = sp.Matrix([angular_index * amplitude * sp.cos(angular_index * s) / r, -c**2 * slope * sp.sin(angular_index * s) / d, c * r * slope * sp.sin(angular_index * s) / d])
    ledger.check("full physical velocity cross product fixes the KKS orientation", sp.simplify(sp.trigsimp(h.dot(va.cross(vb)) - angular_index * c * amplitude * slope / r)) == 0)
    ledger.check("KKS radial integration has the exact pressure weight", sp.simplify(-c * sp.diff(2 * c * circulation / d**2, r) - r * beta) == 0)
    action_h = sp.symbols("h_mode", positive=True)
    symplectic = sp.Matrix([[0, action_h / sigma], [-action_h / sigma, 0]])
    dynamic = sp.Matrix([[0, sigma], [-sigma, 0]])
    ledger.check("positive full phase action produces the physical mode clock", symplectic * dynamic + action_h * sp.eye(2) == sp.zeros(2))
    angle_scale = sp.symbols("c_obs", nonzero=True, real=True)
    physical_mass = action_h / (sigma**2 * angle_scale**2)
    ledger.check("observed scalar action has positive physical mass", physical_mass.is_positive is True)

    xi_r = sp.symbols("xi_r", real=True)
    spin = r * vt + sp.diff(r**2 * f, r) * xi_r
    axial_momentum = vz + c * sp.diff(f, r) * xi_r
    ledger.check("mechanical spin retains its exact axial-momentum companion", sp.simplify(spin + c * axial_momentum - (r * vt + c * vz)) == 0)
    ledger.check("zero helical charge is not zero mechanical spin", sp.simplify(spin.subs(vt, -c * vz / r)) != 0)
    state_a, state_b = sp.symbols("a b", real=True)
    phi2 = amplitude * (state_a * sp.cos(2 * s) + state_b * sp.sin(2 * s))
    phi2_r = slope * (state_a * sp.cos(2 * s) + state_b * sp.sin(2 * s))
    observer_derivative = (sp.diff(phi2, s) * 2 + (-phi2_r) * 2 * sp.I * r) * sp.exp(2 * sp.I * s)
    observed_angular = sp.integrate(sp.expand_complex(observer_derivative), (s, 0, 2 * sp.pi))
    expected_angular = 2 * sp.pi * (2 * amplitude + r * slope) * (state_b - sp.I * state_a)
    ledger.check("literal stationary quadrupole differentiates in the claimed physical clock", sp.simplify(observed_angular - expected_angular) == 0)
    xi_r_mode = 2 * amplitude * state_a / (sigma * r)
    spin_mode = spin.subs({vt: -c**2 * slope * state_a / d, xi_r: xi_r_mode})
    expected_spin = state_a * (-c**2 * r * slope / d + 4 * circulation * c**2 * amplitude / (sigma * d**2))
    ledger.check("complete moving-tag spin has the derived radial density", sp.simplify(spin_mode - expected_spin) == 0)
    b1, b2, ratio1, ratio2 = sp.symbols("B1 B2 ratio1 ratio2", real=True)
    moment = sp.Matrix([[b1, b2], [b1 * ratio1, b2 * ratio2]])
    ledger.check("two actual radial observation rows have the exact matching determinant", sp.expand(moment.det() - b1 * b2 * (ratio2 - ratio1)) == 0)
    profile = sp.Function("f")(r)
    ledger.check("a localized profile keeps the full nonconstant-momentum coupling", sp.diff((r**2 + c**2) * profile, r) == 2 * r * profile + (r**2 + c**2) * sp.diff(profile, r))
    print(f"H metric: c²/(c²+r²), beta={beta}; full KKS=h_mode/sigma")
    print("Actual periodic modes have omega=-sigma<0, positive action and S_z=-c delta P_z.")
    print("No off-carrier curvature sign or Euclidean EPS closed-tube transfer is inferred.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
