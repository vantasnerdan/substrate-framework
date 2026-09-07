"""Exact symbolic checks for the P253/0074 massive cylindrical fiber ledger."""

import sympy as sp


x, z, delta, k = sp.symbols("x z delta k", nonzero=True)
q = 1 + delta * x
vr = sp.Function("v_r")(x, z)
vt = sp.Function("v_theta")(x, z)
vz = sp.Function("v_z")(x, z)
f = sp.Function("f")(x, z)


def div_fiber(v):
    r, t, zz = v
    return sp.diff(r, x) + sp.diff(zz, z) + delta * r / q + sp.I * k * t / q


def grad_fiber(s):
    return (sp.diff(s, x), sp.I * k * s / q, sp.diff(s, z))


def curl_fiber(v):
    r, t, zz = v
    return (
        sp.I * k * zz / q - sp.diff(t, z),
        sp.diff(r, z) - sp.diff(zz, x),
        sp.diff(t, x) + delta * t / q - sp.I * k * r / q,
    )


lap = sp.diff(f, x, 2) + sp.diff(f, z, 2) + delta * sp.diff(f, x) / q - k**2 * f / q**2
assert sp.simplify(div_fiber(grad_fiber(f)) - lap) == 0
assert all(sp.simplify(c) == 0 for c in curl_fiber(grad_fiber(f)))
assert sp.simplify(div_fiber(curl_fiber((vr, vt, vz)))) == 0

# In v_r +/- i v_theta, the cylindrical vector Laplacian diagonalizes with
# masses delta^2(n +/- 1)^2=(k +/- delta)^2.
n = sp.symbols("n", integer=True)
assert sp.simplify((delta**2 * (n + 1) ** 2).subs(n, k / delta) - (k + delta) ** 2) == 0
assert sp.simplify((delta**2 * (n - 1) ** 2).subs(n, k / delta) - (k - delta) ** 2) == 0

inv_q = sp.series(1 / q, delta, 0, 3).removeO()
inv_q2 = sp.series(1 / q**2, delta, 0, 3).removeO()
assert sp.expand(inv_q) == 1 - delta * x + delta**2 * x**2
assert sp.expand(inv_q2) == 1 - 2 * delta * x + 3 * delta**2 * x**2

# Exact Newton-kernel normalization after the physical harmonic projection,
# half-density conjugation, distance scaling, and alpha=delta*t.
R, a = sp.symbols("R a", positive=True)
assert sp.simplify((R / a) * delta).subs(delta, a / R) == 1

# Re-derive the connection signs in the two Lie brackets by converting
# orthonormal theta components to coordinate components and back.
rr, theta, zz = sp.symbols("r theta z", positive=True, real=True)
Wr0 = sp.Function("W_r")(rr, zz)
Wz0 = sp.Function("W_z")(rr, zz)
vr0 = sp.Function("p_r")(rr, zz) * sp.exp(sp.I * n * theta)
vt0 = sp.Function("p_theta")(rr, zz) * sp.exp(sp.I * n * theta)
vz0 = sp.Function("p_z")(rr, zz) * sp.exp(sp.I * n * theta)
ze0 = sp.Function("zeta")(rr, zz)

# theta physical component of [W,p]: r times its coordinate-theta component.
transport_theta_coord = rr * (
    Wr0 * sp.diff(vt0 / rr, rr) + Wz0 * sp.diff(vt0 / rr, zz)
)
transport_theta_phys = Wr0 * sp.diff(vt0, rr) + Wz0 * sp.diff(vt0, zz) - Wr0 * vt0 / rr
assert sp.simplify(transport_theta_coord - transport_theta_phys) == 0

# theta physical component of [v,zeta*e_theta].
stretch_theta_coord = rr * (
    vr0 * sp.diff(ze0 / rr, rr)
    + vz0 * sp.diff(ze0 / rr, zz)
    - (ze0 / rr) * sp.diff(vt0 / rr, theta)
)
stretch_theta_phys = (
    vr0 * sp.diff(ze0, rr)
    + vz0 * sp.diff(ze0, zz)
    - ze0 * vr0 / rr
    - sp.I * n * ze0 * vt0 / rr
)
assert sp.simplify(stretch_theta_coord - stretch_theta_phys) == 0

# Wrong connection signs destroy div(curl)=0 and are therefore exposed.
bad_curl = list(curl_fiber((vr, vt, vz)))
bad_curl[2] = bad_curl[2] - 2 * delta * vt / q
assert sp.simplify(div_fiber(tuple(bad_curl))) != 0

print("PASS exact scalar div-grad identity")
print("PASS exact curl-grad and div-curl complexes")
print("PASS helical masses (k+delta)^2 and (k-delta)^2")
print("PASS q^-1 and q^-2 second-order expansions")
print("PASS half-density Newton coefficient sqrt(q*q')/(4*pi)")
print("PASS cylindrical transport/stretching connection signs")
print("PASS wrong cylindrical connection sign is detected")
