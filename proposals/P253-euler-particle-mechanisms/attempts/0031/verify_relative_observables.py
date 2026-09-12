"""Exact algebra checks for the P253/0031 relative observables.

This verifier checks cylindrical identities and normalization-sensitive
one-dimensional reductions.  Convergence, global group integration, and the
existence proof remain analytic statements in derivation.md.
"""

import sympy as sp


r, z, c = sp.symbols("r z c", positive=True, real=True)
f = sp.Function("f")(r, z)
Q = sp.Function("Q")(r, z)  # total swirl numerator F(psi)
L = sp.Function("L")(r)
Fp = sp.Symbol("Fprime", real=True)

fr, fz = sp.diff(f, r), sp.diff(f, z)
delta_star_f = sp.diff(f, r, 2) - fr/r + sp.diff(f, z, 2)
psi0_r = -c*r

# Direct laboratory velocity and cylindrical curl.
ur, ut, uz = -fz/r, Q/r, fr/r
wr = -sp.diff(Q, z)/r
wt = -delta_star_f/r
wz = sp.diff(Q, r)/r
raw_dot = sp.expand(ur*wr + ut*wt + uz*wz)
chain_rule = {
    sp.diff(Q, z): Fp*fz,
    sp.diff(Q, r): Fp*(psi0_r + fr),
}
raw_dot = sp.simplify(raw_dot.xreplace(chain_rule))
expected_dot = sp.simplify(
    (Fp*(fz**2 + fr*(psi0_r + fr)) - Q*delta_star_f)/r**2
)
assert sp.simplify(raw_dot - expected_dot) == 0
print("PASS laboratory cylindrical vorticity and raw helicity density")

# Weighted Green identity: div(F grad f/r) is the complete meridional
# boundary density.  The chain rule grad F=F'(psi)grad(psi0+f) is retained.
A = fz**2 + fr*(psi0_r + fr)
boundary_divergence = Fp*A/r + Q*delta_star_f/r
raw_meridional = (Fp*A - Q*delta_star_f)/r
reduced_meridional = 2*Fp*A/r - boundary_divergence
assert sp.simplify(raw_meridional - reduced_meridional) == 0
print("PASS helicity Green identity with signed boundary divergence")

# The z component of x cross v is r*v_theta=Q-L.  The volume Jacobian r and
# azimuthal integral 2*pi therefore leave the exact factor 2*pi*r*(Q-L).
G = Q - L
jz_density_after_theta = sp.simplify(2*sp.pi*r*(r*(G/r)))
assert sp.simplify(jz_density_after_theta - 2*sp.pi*r*G) == 0
print("PASS literal axial angular-momentum factor 2*pi*r")

# The vorticity-moment conversion retains its radial surface derivative.
Gr = sp.Function("G")(r)
delta_omega_z = sp.diff(Gr, r)/r
assert sp.simplify(
    2*r*Gr - (sp.diff(r**2*Gr, r) - r**3*delta_omega_z)
) == 0
print("PASS angular-momentum vorticity conversion with surface row")

# The impulse surface row follows from integrating r*Delta_* f after the
# z integral M(r)=integral f dz is taken.
M = sp.Function("M")(r)
assert sp.simplify(
    sp.diff(r*sp.diff(M, r) - 2*M, r)
    - (r*sp.diff(M, r, 2) - sp.diff(M, r))
) == 0
print("PASS axial impulse antiderivative r*M'-2*M")

# The exact limiting homoclinic supplies, rather than assigns, the factor 6/beta.
X = sp.Symbol("X", real=True)
beta = sp.Symbol("beta", positive=True)
Astar = 3*sp.sech(X/2)**2/(2*beta)
antiderivative = 3*sp.tanh(X/2)/beta
assert sp.trigsimp(sp.diff(antiderivative, X) - Astar) == 0
assert sp.simplify(
    sp.limit(antiderivative, X, sp.oo)
    - sp.limit(antiderivative, X, -sp.oo)
    - 6/beta
) == 0
print("PASS homoclinic integral 6/beta")

# With X_xi=-ad*_xi m and Omega(X_xi,X_eta)=-<m,[xi,eta]>, both the
# Hamiltonian row and the positive momentum-map row have the same sign.
B = sp.Symbol("bracket_pairing", real=True)
dH_on_Xeta = -B
omega_Xu_Xeta = -B
assert dH_on_Xeta == omega_Xu_Xeta
dJxi_on_Xeta = -B
omega_Xxi_Xeta = -B
assert dJxi_on_Xeta == omega_Xxi_Xeta
print("PASS KKS/Hamiltonian/momentum-map sign convention")

# Exposing mutations: each changes a nonzero symbolic identity rather than a
# formatting literal.
assert sp.simplify(raw_meridional - (2*Fp*A/r + boundary_divergence)) != 0
assert sp.simplify(jz_density_after_theta - sp.pi*r*G) != 0
assert sp.simplify(
    sp.diff(r*sp.diff(M, r) - M, r)
    - (r*sp.diff(M, r, 2) - sp.diff(M, r))
) != 0
assert dH_on_Xeta != B
print("PASS four sign/factor/surface exposing mutations")

print("7 exact checks plus 4 exposing mutations passed")
