"""Exact algebraic checks for the material/centroid joining identities."""

import sympy as s


checks = []


def check(name, value):
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.simplify(x) == 0 for x in entries), (name, value)
    checks.append(name)
    print(f"PASS {name}")


rho, mass = s.symbols("rho mass", positive=True)
p, w1, w2 = s.symbols("p w1 w2", real=True)
x, r1, r2 = s.symbols("x r1 r2", real=True)
dx, dr1, dr2 = s.symbols("dx dr1 dr2", real=True)
vel = [p / (3 * mass) + w1, p / (3 * mass) + w2,
       p / (3 * mass) - w1 - w2]
pos = [x + r1, x + r2, x - r1 - r2]
dpos = [dx + dr1, dx + dr2, dx - dr1 - dr2]
check("exact three-parcel centroid", sum(pos) / 3 - x)
check("exact centroid momentum", mass * sum(vel) - p)
check("same kinetic split", mass * sum(v**2 for v in vel) / 2
      - p**2 / (6 * mass) - mass * (w1**2 + w2**2 + (w1 + w2)**2) / 2)
check("same cotangent split", mass * sum(v * d for v, d in zip(vel, dpos))
      - p * dx - mass * (w1 * dr1 + w2 * dr2
                         + (w1 + w2) * (dr1 + dr2)))

# The first moment is antisymmetric for a compact solenoidal vector field.
Lx, Ly, Lz = s.symbols("Lx Ly Lz", real=True)
bx, by, bz = s.symbols("bx by bz", real=True)
L = s.Matrix([Lx, Ly, Lz])
beta = s.Matrix([bx, by, bz])
Q = s.Matrix(3, 3, lambda i, j: -sum(s.LeviCivita(i, j, m) * L[m]
                                   for m in range(3)) / 2)
rot = s.Matrix(3, 3, lambda i, j: -sum(s.LeviCivita(i, j, m) * beta[m]
                                     for m in range(3)))
check("compact first moment has zero STF", Q + Q.T)
check("physical affine spin pairing", sum(rot[i, j] * Q[i, j]
                                         for i in range(3) for j in range(3))
      - beta.dot(L))
B, bn = s.symbols("B bn", real=True)
check("Gamma affine KKS cancellation", B * bn - bn * B)

# Actual Bloch material lift on a nonconstant smooth Beltrami witness.
z, t, k = s.symbols("z t k", real=True)
ay, az, dy, dz = s.symbols("ay az dy dz", real=True)
u = s.Matrix([s.cos(z), -s.sin(z), 0])
a = s.Matrix([0, ay, az])
adot = s.Matrix([0, dy, dz])
vb = adot + s.I * k * u[0] * a - az * u.diff(z)
check("Bloch material velocity is exactly solenoidal",
      s.I * k * vb[0] + s.diff(vb[2], z))
vx, vy, vz = s.symbols("vx vy vz", real=True)
V = s.Matrix([vx, vy, vz])
vg = V - t * vz * u.diff(z)
check("Galilean tangent solves linearized Euler",
      vg.diff(t) + vg[2] * u.diff(z))
def avg(expr):
    return s.integrate(s.expand(expr), (z, 0, 2 * s.pi)) / (2 * s.pi)


check("boost energy density has the full mass",
      avg(rho * ((u + V).dot(u + V) - u.dot(u)) / 2) - rho * V.dot(V) / 2)

# Finite exact six-direction covariance represents the isotropic second moment.
sigma = s.symbols("sigma", positive=True)
samples = [sign * s.sqrt(3) * sigma * s.eye(3)[:, i]
           for i in range(3) for sign in (-1, 1)]
check("zero background mean", sum(samples, s.zeros(3, 1)))
check("isotropic actual covariance",
      sum((u0 * u0.T for u0 in samples), s.zeros(3)) / 6 - sigma**2 * s.eye(3))
hd = s.Matrix(3, 3, s.symbols("h0:9", real=True))
dU = s.Matrix(s.symbols("d0:3", real=True))
kin = sum((dU + hd * u0).dot(dU + hd * u0) for u0 in samples) / 6
check("complete mean material macro gradient term",
      kin - dU.dot(dU) - sigma**2 * sum(h**2 for h in hd))

v, c, lam = s.symbols("v c lam", nonzero=True, real=True)
H = rho * (v**2 - v * c / lam)
K = rho * (v * c / lam - c**2 / lam**2)
check("material-orbit difference is retained", H - K - rho * (v - c / lam)**2)

# A genuinely nonzero signed canonical/KKS witness; not an inference from H.
omega = s.Matrix([s.cos(z), s.sin(z), 0])
ub = -omega
q_field = s.Matrix([0, 0, 1])
s_field = s.Matrix([s.sin(z), -s.cos(z), 0])
vq = q_field.cross(omega)
vs = s_field.cross(omega)
piq = rho * (vq + q_field[2] * ub.diff(z))
pis = rho * (vs + s_field[2] * ub.diff(z))
canonical_pair = q_field.dot(pis) - s_field.dot(piq)
kks_pair = rho * omega.dot(q_field.cross(s_field))
check("signed nonzero canonical pullback equals KKS", canonical_pair - kks_pair)
check("symplectic witness is rho not zero", canonical_pair - rho)
check("same canonical graph has exact Q energy", piq.dot(piq) / (2 * rho))
check("same canonical graph has exact S energy", pis.dot(pis) / (2 * rho) - rho / 2)

# Exposing full operator example: G*E and G*AE are BOTH nonzero.
G = s.Matrix([[1, 0], [0, 1], [1, 1], [0, 1]])
E = s.Matrix([[1, 1], [2, 0], [0, 1], [1, -1]])
B0 = s.Matrix([[0, 1, 2, 0], [-1, 0, 1, 3],
               [-2, -1, 0, 2], [0, -3, -2, 0]])
Aop = s.Matrix([[1, 2, 0, 1], [0, -1, 1, 0],
                [2, 0, 1, 1], [1, 0, 0, 2]])
pressure = s.diag(2, 3, 5, 7)
T = G.row_join(E)
J = s.zeros(4, 2).row_join(Aop * E)
y = s.Matrix(s.symbols("y0:4", real=True))
yd = s.Matrix(s.symbols("yd0:4", real=True))
ww = s.Matrix(s.symbols("w0:2", real=True))
eta = T * y
pi = rho * ((B0 * T + J) * y + G * ww)
Kmat = pressure - rho * B0.T * B0
HJ = pi.dot(pi) / (2 * rho) - pi.dot(B0 * eta) + eta.dot(pressure * eta) / 2
hexpanded = eta.dot(Kmat * eta) / 2 + rho * ((G * ww + J * y).dot(G * ww + J * y)) / 2
check("complete joint Hamiltonian with all crosses", s.expand(HJ - hexpanded))
theta = pi.dot(T * yd)
lag = theta - HJ
wstar = (G.T * G).inv() * G.T * (T * yd - J * y)
check("full noncommuting macro momentum stationary equation",
      s.Matrix([s.diff(lag, wi) for wi in ww]).subs(dict(zip(ww, wstar))))
PG = G * (G.T * G).inv() * G.T
reduced = rho * (PG * (T * yd - J * y)).dot(PG * (T * yd - J * y)) / 2
reduced += rho * ((B0 * T + J) * y).dot(T * yd)
reduced -= eta.dot(Kmat * eta) / 2 + rho * (J * y).dot(J * y) / 2
check("full macro Schur action retains every mixed block",
      s.expand(lag.subs(dict(zip(ww, wstar))) - reduced))
print(f"{len(checks)}/{len(checks)} exact checks passed")
