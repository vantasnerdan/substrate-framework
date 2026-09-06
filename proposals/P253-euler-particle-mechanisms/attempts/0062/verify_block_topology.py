#!/usr/bin/env python3
"""Exact algebra for the 0062 DA block, centralizer, and layer ledger.

The checks prove the displayed local identities only.  They do not construct
the missing global Cao graph domain, limiting-absorption resolvent, mode, or
nonlinear branch.
"""

import sympy as sp


r, theta, z = sp.symbols("r theta z", positive=True, real=True)
n = sp.symbols("n", integer=True, nonzero=True)
Omega = sp.symbols("Omega", real=True)
phase = sp.exp(sp.I * n * theta)


def bracket(X, Y):
    """Coordinate-basis Lie bracket in (r,theta,z)."""
    coords = (r, theta, z)
    return sp.Matrix(
        [
            sp.simplify(
                sum(X[j] * sp.diff(Y[i], coords[j]) for j in range(3))
                - sum(Y[j] * sp.diff(X[i], coords[j]) for j in range(3))
            )
            for i in range(3)
        ]
    )


Wr = sp.Function("W_r")(r, z)
Wz = sp.Function("W_z")(r, z)
zeta = sp.Function("zeta")(r, z)
W = sp.Matrix([Wr, -Omega, Wz])
omega = sp.Matrix([0, zeta, 0])

er = sp.Function("eta_r")(r, z)
et = sp.Function("eta_theta")(r, z)
ez = sp.Function("eta_z")(r, z)
hr = sp.Function("h_r")(r, z)
ht = sp.Function("h_theta")(r, z)
hz = sp.Function("h_z")(r, z)
eta = phase * sp.Matrix([er, et, ez])
hvel = phase * sp.Matrix([hr, ht, hz])

T_er = Wr * sp.diff(er, r) + Wz * sp.diff(er, z) - sp.I * n * Omega * er
T_et = Wr * sp.diff(et, r) + Wz * sp.diff(et, z) - sp.I * n * Omega * et
T_ez = Wr * sp.diff(ez, r) + Wz * sp.diff(ez, z) - sp.I * n * Omega * ez
expected_A = phase * sp.Matrix(
    [
        -T_er + er * sp.diff(Wr, r) + ez * sp.diff(Wr, z) + sp.I * n * zeta * hr,
        -T_et - hr * sp.diff(zeta, r) - hz * sp.diff(zeta, z) + sp.I * n * zeta * ht,
        -T_ez + er * sp.diff(Wz, r) + ez * sp.diff(Wz, z) + sp.I * n * zeta * hz,
    ]
)
actual_A = -bracket(W, eta) - bracket(hvel, omega)
assert sp.simplify(actual_A - expected_A) == sp.zeros(3, 1)

# Exact nonzero-harmonic DA map and its positive-core inverse.
xr, xt, xz = sp.symbols("xi_r xi_theta xi_z")
xi = phase * sp.Matrix([xr, xt, xz])
Cxi = sp.simplify(-bracket(xi, omega) / phase)
expected_C = sp.Matrix(
    [
        sp.I * n * zeta * xr,
        -xr * sp.diff(zeta, r) - xz * sp.diff(zeta, z) + sp.I * n * zeta * xt,
        sp.I * n * zeta * xz,
    ]
)
assert sp.simplify(Cxi - expected_C) == sp.zeros(3, 1)

ar, at, az = sp.symbols("a_r a_theta a_z")
xi_r_inv = ar / (sp.I * n * zeta)
xi_z_inv = az / (sp.I * n * zeta)
xi_t_inv = (at + xi_r_inv * sp.diff(zeta, r) + xi_z_inv * sp.diff(zeta, z)) / (
    sp.I * n * zeta
)
recovered = expected_C.subs({xr: xi_r_inv, xt: xi_t_inv, xz: xi_z_inv})
assert sp.simplify(recovered - sp.Matrix([ar, at, az])) == sp.zeros(3, 1)

# Axisymmetric positive-core centralizer: psi=F(zeta) gives a divergence-free
# poloidal field tangent to every regular zeta level.
F = sp.Function("F")
psi = F(zeta)
Yr = -sp.diff(psi, z) / r
Yz = sp.diff(psi, r) / r
tangency = sp.simplify(Yr * sp.diff(zeta, r) + Yz * sp.diff(zeta, z))
divergence = sp.simplify(sp.diff(r * Yr, r) / r + sp.diff(Yz, z))
assert tangency == 0
assert divergence == 0

# Physical KKS density. omega=zeta*partial_theta=r*zeta*e_theta and
# dV=r dr dtheta dz, hence the poloidal density is r^2*zeta.
xrp, xtp, xzp, yrp, ytp, yzp = sp.symbols("x_r x_t x_z y_r y_t y_z")
omega_phys = sp.Matrix([0, r * zeta, 0])
xi_phys = sp.Matrix([xrp, r * xtp, xzp])
chi_phys = sp.Matrix([yrp, r * ytp, yzp])
kks_integrand = sp.expand(r * omega_phys.dot(xi_phys.cross(chi_phys)))
expected_kks = r**2 * zeta * (xzp * yrp - xrp * yzp)
assert sp.simplify(kks_integrand - expected_kks) == 0

# An exact divergence-free inner generator for volume form J dI dphi dsigma.
I, phi, sigma = sp.symbols("I phi sigma", real=True)
J = sp.Function("J")(I, phi, sigma)
S = sp.Function("S")(I, phi, sigma)
XI = sp.diff(S, phi) / J
Xphi = -sp.diff(S, I) / J
volume_divergence = sp.simplify(
    (sp.diff(J * XI, I) + sp.diff(J * Xphi, phi)) / J
)
assert volume_divergence == 0

# Cutoff scaling exponents for layer width h and amplitude h^a.
a, order, sobolev = sp.symbols("a order s", real=True)
l2_derivative_exponent = sp.simplify(a + sp.Rational(1, 2) - order)
assert l2_derivative_exponent.subs({a: 1, order: sobolev}) == sp.Rational(3, 2) - sobolev
assert l2_derivative_exponent.subs({a: 1, order: sobolev + 1}) == sp.Rational(1, 2) - sobolev
assert l2_derivative_exponent.subs({a: 1, order: sobolev - 1}) == sp.Rational(5, 2) - sobolev

# A genuine pendulum normal form has center/saddle determinants with opposite
# signs, hence a separatrix when alpha*epsilon*V is nonzero.
alpha, epsilon, V = sp.symbols("alpha epsilon V", nonzero=True, real=True)
x = sp.symbols("x", real=True)
H = alpha * x**2 / 2 + epsilon * V * sp.cos(phi)
hessian = sp.hessian(H, (x, phi))
det_zero = sp.simplify(hessian.det().subs(phi, 0))
det_pi = sp.simplify(hessian.det().subs(phi, sp.pi))
assert det_zero == -V * alpha * epsilon
assert det_pi == V * alpha * epsilon
assert sp.simplify(det_zero + det_pi) == 0

print("harmonic_full_hodge_coordinate_block: PASS")
print("positive_core_DA_map_and_inverse: PASS")
print("centralizer_tangency:", tangency)
print("centralizer_divergence:", divergence)
print("KKS_density:", expected_kks)
print("exact_volume_generator_divergence:", volume_divergence)
print("cutoff_Hs_exponent_Aeqh:", sp.Rational(3, 2) - sobolev)
print("cutoff_Hsplus1_exponent_Aeqh:", sp.Rational(1, 2) - sobolev)
print("induced_vorticity_Hsminus1_exponent_Aeqh:", sp.Rational(5, 2) - sobolev)
print("pendulum_hessian_determinants:", det_zero, det_pi)
print("checks: 13 passed")
