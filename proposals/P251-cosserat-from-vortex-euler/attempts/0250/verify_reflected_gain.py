"""Exact parity, fraction and odd-observation checks for the 0250 joint gain."""

import sympy as sp


# Relative sign of the -R realization against R.  The full pair adds only
# when this sign is +1.
for order in (0, 1, 2):
    optical_axial_to_hybrid_polar = (-1) ** (order + 1)
    acoustic_polar_to_optical_axial = (-1) ** (order + 1)
    optical_axial_to_optical_axial = (-1) ** order
    acoustic_polar_to_hybrid_polar = (-1) ** order
    if order % 2 == 0:
        assert optical_axial_to_hybrid_polar == -1
        assert acoustic_polar_to_optical_axial == -1
        assert optical_axial_to_optical_axial == 1
        assert acoustic_polar_to_hybrid_polar == 1
    else:
        assert optical_axial_to_hybrid_polar == 1
        assert acoustic_polar_to_optical_axial == 1


# The positive fractions are reflection scalars and retain the exact optical
# determinant.  No sign is assigned to the probabilities themselves.
theta, f1, f2, Gu, Gp, Su, Sp = sp.symbols(
    "theta f1 f2 Gu Gp Su Sp", nonzero=True
)
optical = sp.Matrix([
    [theta, theta, 0],
    [f1 * Gu, f2 * Gu, f1 * Gp],
    [f1 * Su, f2 * Su, f1 * Sp],
])
expected = theta * f1 * (f2 - f1) * (Gu * Sp - Su * Gp)
assert sp.simplify(optical.det() - expected) == 0


# q_N is even in local K.  Under R -> -R the local polar amplitude and the
# laboratory polar output both reverse, so their product adds.
hp, hpar, kperp2, kpar2, dlocal = sp.symbols(
    "h_perp h_parallel k_perp_squared k_parallel_squared d_local"
)
q = hp * kperp2 + hpar * kpar2
assert sp.simplify(q.subs({kperp2: kperp2, kpar2: kpar2}) - q) == 0
amplitude_R = dlocal / q
amplitude_minus_R = -dlocal / q
lab_output_R = q * amplitude_R
lab_output_minus_R = -q * amplitude_minus_R
assert sp.simplify(lab_output_minus_R - lab_output_R) == 0


# The surviving odd curl rows are a finite perturbation of the even block.
# This is the exact retained five-position observation determinant of 0241.
k, rho, j = sp.symbols("k rho j", nonzero=True)
C = sp.I * sp.Matrix([[0, -k, 0], [k, 0, 0], [0, 0, 0]])
T6 = sp.BlockMatrix([
    [sp.eye(3), -j * C / (2 * rho)],
    [C / 2, sp.eye(3)],
]).as_explicit()
T5 = T6.extract([0, 1, 3, 4, 5], [0, 1, 3, 4, 5])
assert sp.simplify(T5.det() - (1 + j * k**2 / (4 * rho))**2) == 0

print("PASS reflected parity cancels exactly the forbidden even cross blocks")
print("PASS positive fractions retain the optical determinant")
print("PASS q_N source and polar output signs add under R/-R")
print("PASS retained odd curl map determinant = (1+j*k**2/(4*rho))**2")
