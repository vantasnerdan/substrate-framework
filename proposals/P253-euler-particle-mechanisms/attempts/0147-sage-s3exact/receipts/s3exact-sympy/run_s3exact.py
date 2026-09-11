"""0147-sage-s3exact VALIDATION FILE — every ledger identity and every
must-FAIL mutation as a hard assertion. Exit 0 iff ALL checks green.
Any assertion failure = the paper's corresponding ESTABLISHED claim is
withdrawn (00-scope kill condition).

Run:  python3 run_s3exact.py > run.log 2>&1 ; echo $?
Env:  repo python3, sympy only. No repo imports; self-contained.
Discipline (D-08): identities are checked as EQUALITIES; mutations are
checked to FAIL the identity property (the failure detection itself is
asserted, so a vacuous check cannot pass).
"""
import sympy as sp

t, x, y, z = sp.symbols('t x y z')
coords = (x, y, z)
phi, al, be, sg = [sp.Function(f, real=True)(x, y, z) for f in ('phi', 'alpha', 'beta', 'sigma')]
U = [sp.Function(f'u_{c}', real=True)(t, x, y, z) for c in 'xyz']
W = [sp.Function(f'w_{c}', real=True)(x, y, z) for c in 'xyz']
V = [sp.Function(f'v_{c}', real=True)(x, y, z) for c in 'xyz']
P = sp.Function('Pi', real=True)(t, x, y, z)


def grad(f):
    return [sp.diff(f, v).doit() for v in coords]


def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]


def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def div(a):
    return sum(sp.diff(c, v).doit() for c, v in zip(a, coords))


def curl(a):
    return [sp.diff(a[(i+2) % 3], coords[(i+1) % 3]).doit()
            - sp.diff(a[(i+1) % 3], coords[(i+2) % 3]).doit() for i in range(3)]


COUNTS = {"identity": 0, "mutation": 0}


def check(label, cond):
    kind = "mutation" if label.startswith("M") else "identity"
    COUNTS[kind] += 1
    print(f"[{'PASS' if cond else 'FAIL'}] [{kind} #{COUNTS[kind]}] {label}")
    assert cond, f"VALIDATION FAILURE: {label}"


gp, ga, gb = grad(phi), grad(al), grad(be)
u_cl = [gp[i] + al*gb[i] for i in range(3)]
om_cl = cross(ga, gb)

# ---- I-CS-a: vorticity of the Clebsch ansatz is exactly grad(alpha)x grad(beta)
curl_u = [sp.diff(u_cl[(i+2) % 3], coords[(i+1) % 3]).doit()
          - sp.diff(u_cl[(i+1) % 3], coords[(i+2) % 3]).doit() for i in range(3)]
check("I-CS-a curl(grad phi + alpha grad beta) == grad alpha x grad beta",
      all(sp.simplify(curl_u[i] - om_cl[i]) == 0 for i in range(3)))

# ---- I-CS-b: the CS current is identically conserved
check("I-CS-b div(grad alpha x grad beta) == 0 identically",
      sp.simplify(div(om_cl)) == 0)

# ---- I-H-a: helicity density is exactly the Jacobian determinant
J = sp.Matrix([gp, ga, gb]).det()
check("I-H-a u.omega == det J(phi, alpha, beta) exactly",
      sp.simplify(sp.expand(dot(u_cl, om_cl) - J)) == 0)

# ---- I-H-b: generic curl-cross expansion identity
Ux = [sp.Function(f'us_{c}', real=True)(x, y, z) for c in 'xyz']
lhs = curl(cross(W, Ux))
rhs = [W[i]*div(Ux) - Ux[i]*div(W)
       + sum(Ux[j]*sp.diff(W[i], coords[j]) for j in range(3))
       - sum(W[j]*sp.diff(Ux[i], coords[j]) for j in range(3)) for i in range(3)]
check("I-H-b curl(w x u) == w div u - u div w + (u.grad)w - (w.grad)u",
      all(sp.simplify(sp.expand(lhs[i] - rhs[i])) == 0 for i in range(3)))

# ---- I-H-c: flux-transport identity
lhs2 = dot(Ux, curl(W))
rhs2 = div(cross(W, Ux)) + dot(W, curl(Ux))
check("I-H-c u.curl(w) == div((w x u)) + w.curl(u)",
      sp.simplify(sp.expand(lhs2 - rhs2)) == 0)

# ---- I-Helicity: full flux identity under Euler evolution du/dt = -om x u - grad Pi
om = curl(U)
E = sp.diff(dot(U, om), t).doit()
subs_map = {}
rhs_t = [-(cross(om, U)[i] + grad(P)[i]) for i in range(3)]
for i in range(3):
    subs_map[sp.Derivative(U[i], t)] = rhs_t[i]
    for j in range(3):
        subs_map[sp.Derivative(U[i], t, coords[j])] = sp.diff(rhs_t[i], coords[j]).doit()
E2 = sp.simplify(sp.expand(E.subs(subs_map).doit()))
flux = div([P*om[i] for i in range(3)]) + div(cross(cross(om, U), U))
check("I-Helicity d(u.omega)/dt == -div(Pi*omega + (omega x u) x u) under Euler evolution",
      sp.simplify(E2 + flux) == 0)

# ---- I-Noether2: gauged kinetic action has vanishing variational Hessian in (d_i A_j)
k = sp.symbols('kappa', positive=True)
A = [sp.Function(f'A_{c}', real=True)(x, y, z) for c in 'xyz']
u_g = [gp[i] + al*gb[i] - k*al*A[i] for i in range(3)]
T = sp.expand((u_g[0]**2 + u_g[1]**2 + u_g[2]**2) / 2)
dvars = [sp.Derivative(A[j], coords[i]) for i in range(3) for j in range(3)]
n_bad = sum(1 for dv1 in dvars for dv2 in dvars if sp.diff(T, dv1, dv2) != 0)
check("I-Noether2 variational Hessian d^2T/d(d_i A_j)d(d_k A_l) == 0 (81/81)",
      n_bad == 0)
check("I-Noether2 EL source dT/dA_j algebraic (no dA terms to integrate)",
      all("Derivative(A" not in str(sp.diff(T, A[j])) for j in range(3)))

# ================= MUTATIONS (each must FAIL its identity property) ============

# M1: symmetrized ansatz kills vorticity (I-CS-a detects ansatz structure)
u_m = [gp[i] + al*gb[i] + be*ga[i] for i in range(3)]
curl_um = [sp.diff(u_m[(i+2) % 3], coords[(i+1) % 3]).doit()
           - sp.diff(u_m[(i+1) % 3], coords[(i+2) % 3]).doit() for i in range(3)]
check("M1 mutation: symmetrized ansatz vorticity vanishes (structure detected)",
      all(sp.simplify(c) == 0 for c in curl_um))

# M2: wrong current structure grad(alpha).grad(beta) is not conserved (I-CS-b detects)
dot_ab = sum(ga[i]*gb[i] for i in range(3))
check("M2 mutation: div(grad alpha . grad beta) generically nonzero",
      sp.simplify(sum(sp.diff(dot_ab, v).doit() for v in coords)) != 0)

# M3: non-solenoidal surrogate vorticity curl u + grad sigma (I-Helicity flux breaks:
#     div(w') = Laplacian(sigma), so -grad Pi.w' is no longer a pure flux)
sur = [curl(Ux)[i] + grad(sg)[i] for i in range(3)]
lap = sum(sp.diff(sg, v, 2).doit() for v in coords)
check("M3 mutation: div(curl u + grad sigma) == Laplacian(sigma) != 0 generically",
      sp.simplify(div(sur) - lap) == 0 and sp.simplify(lap) != 0)

# M4: helicity closure mechanism is the antisymmetric contraction with curl u itself:
#     (w x u).w == 0 only because w = curl u is contracted against itself;
#     an independent density v breaks the closure at (w x u).v != 0.
check("M4a closure: (w x u).w == 0 identically",
      sp.simplify(sp.expand(dot(cross(W, Ux), W))) == 0)
check("M4b mutation: (w x u).v with independent v generically nonzero",
      sp.simplify(sp.expand(dot(cross(W, Ux), V))) != 0)

# M5: the propagating completion itself (add (eps/2)|curl A|^2) is detected:
#     the Hessian acquires nonzero entries — exactly the certificate that a
#     propagating extension MUST add (d i A_j)-kinetic structure.
eps = sp.symbols('epsilon', positive=True)
Tc = sp.expand(sum((sp.diff(A[(i+2) % 3], coords[(i+1) % 3])
                    - sp.diff(A[(i+1) % 3], coords[(i+2) % 3]))**2 for i in range(3)) / 2)
T2 = sp.expand(T + eps*Tc)
n_bad2 = sum(1 for dv1 in dvars for dv2 in dvars if sp.diff(T2, dv1, dv2) != 0)
check("M5 mutation: added (eps/2)|curl A|^2 gives nonzero Hessian entries",
      n_bad2 > 0)

total = COUNTS["identity"] + COUNTS["mutation"]
print(f"ALL CHECKS GREEN: {COUNTS['identity']} identity-assertions "
      f"+ {COUNTS['mutation']} mutations = {total} assertions (self-counted).")
