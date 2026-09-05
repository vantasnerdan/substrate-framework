"""Exact physical two-branch pullback against the canonical bulk operator.

K=(0,0,k) is without loss for the declared isotropic tensors. Restrict
only U_z=0; retain Phi_z. rho,j are physical mass/spin density, nu>0.
No Euler supplier is encoded as an assumed numerical matrix equality.
"""

import sympy as s

from substrate_framework.micropolar import (
    MicropolarCoefficients,
    micropolar_fourier_stiffness,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0227-physical-branch-reconstruction")
    rho, j, nu, a, ct, cl = s.symbols("rho j nu a cT cL", positive=True)
    k = s.symbols("k", real=True)
    curl = s.I*s.Matrix([[0, -k, 0], [k, 0, 0], [0, 0, 0]])
    beta = j/(2*rho)
    full = s.BlockMatrix([[s.eye(3), -beta*curl], [curl/2, s.eye(3)]]).as_explicit()
    indices = [0, 1, 3, 4, 5]
    transform = full.extract(indices, indices)
    mass0 = s.diag(rho, rho, j, j, j)
    branch_mass = s.simplify(transform.conjugate().T*mass0*transform)
    expected_mass = s.diag(rho+j*k**2/4, rho+j*k**2/4,
                           j+j**2*k**2/(4*rho), j+j**2*k**2/(4*rho), j)
    checks.check("physical polar/axial observations cancel the complete kinetic cross blocks",
                 branch_mass == expected_mass)
    checks.check("the exact physical map is nonsingular at every real wave number",
                 s.factor(transform.det()-(1+j*k**2/(4*rho))**2) == 0)
    inverse = s.simplify(transform.inv())
    checks.check("both second-gradient branch masses return exactly the physical metric",
                 s.simplify(inverse.conjugate().T*branch_mass*inverse-mass0) == s.zeros(5))
    frequencies = s.diag(a*k**2, a*k**2, nu**2+ct*k**2,
                         nu**2+ct*k**2, nu**2+cl*k**2)
    stiffness = s.simplify(inverse.conjugate().T*branch_mass*frequencies*inverse)
    jet = stiffness.applyfunc(lambda value: s.series(value, k, 0, 3).removeO().expand())
    alpha = j*nu**2/4
    gt, gl = j*(ct-alpha/rho), j*cl
    cs = s.symbols("c_s", real=True)
    coefficients = MicropolarCoefficients(0, rho*a, alpha, gl/2-cs, cs, gt-cs)
    canonical = s.Matrix(micropolar_fourier_stiffness([0, 0, k], coefficients))
    canonical = canonical.extract(indices, indices)
    checks.check("the complete physical second jet equals the canonical coupled operator",
                 s.simplify(jet-canonical) == s.zeros(5))
    checks.check("both physical off-diagonal blocks have the same Hermitian coupling",
                 s.simplify(jet[:2, 2:]+2*alpha*curl[:2, :]) == s.zeros(2, 3)
                 and jet[2:, :2] == jet[:2, 2:].conjugate().T)
    checks.check("the transverse optical curvature retains its translational reaction shift",
                 s.simplify(jet[2, 2].coeff(k, 2)-j*ct+j*alpha/rho) == 0)
    checks.check("incompressible displacement does not remove longitudinal spin",
                 s.simplify(jet[4, 4]-j*(nu**2+cl*k**2)) == 0)
    checks.check("the exact model retains higher spatial terms beyond its Cosserat jet",
                 s.simplify(stiffness-jet) != s.zeros(5))
    checks.check("the exact pulled-back branch pencil has its defining eigenvalues",
                 s.simplify(stiffness*transform-mass0*transform*frequencies) == s.zeros(5))
    checks.check("dropping the second-gradient branch mass changes physical inertia",
                 s.simplify(inverse.conjugate().T*mass0*inverse-mass0) != s.zeros(5))
    wrong = full.subs(beta, 2*beta) if beta.is_Symbol else s.BlockMatrix(
        [[s.eye(3), -2*beta*curl], [curl/2, s.eye(3)]]).as_explicit()
    wrong = wrong.extract(indices, indices)
    wrong_cross = (wrong.conjugate().T*mass0*wrong)[:2, 2:]
    checks.check("doubling the optical current factor leaves a nonzero kinetic cross defect",
                 s.simplify(wrong_cross) != s.zeros(2, 3))
    checks.check("the bulk Fourier operator cannot fix the boundary-divergence modulus",
                 canonical.diff(cs) == s.zeros(5))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
