"""Full Euler residual and bounded-helicity periodic construction identities."""

import sympy as s

from substrate_framework import euler_fourier as fourier
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0145-global-toroidal")
    r, z = s.symbols("r z", real=True, positive=True)
    kap, k, lam, amplitude = s.symbols("kappa k lambda A", positive=True)
    j0, j1 = s.besselj(0, kap*r), s.besselj(1, kap*r)
    psi = amplitude*r*j1*s.cos(k*z)
    ur = amplitude*k*j1*s.sin(k*z)
    ut = amplitude*lam*j1*s.cos(k*z)
    uz = amplitude*kap*j0*s.cos(k*z)
    field = (ur, ut, uz)

    def reduced(value):
        return s.simplify(s.expand(s.expand_func(value)).subs(lam**2, kap**2+k**2))

    ledger.check("explicit global field is the complete streamfunction construction",
                 all(reduced(v) == 0 for v in (ur+s.diff(psi, z)/r,
                                               ut-lam*psi/r, uz-s.diff(psi, r)/r)))
    ledger.check("global cylindrical divergence vanishes without a wall",
                 reduced(s.diff(r*ur, r)/r+s.diff(uz, z)) == 0)
    curl = (-s.diff(ut, z), s.diff(ur, z)-s.diff(uz, r), s.diff(r*ut, r)/r)
    ledger.check("all three curl components have the same exact eigenvalue",
                 all(reduced(curl[i]-lam*field[i]) == 0 for i in range(3)))
    p = -sum(v*v for v in field)/2
    euler = (ur*s.diff(ur, r)+uz*s.diff(ur, z)-ut**2/r+s.diff(p, r),
             ur*s.diff(ut, r)+uz*s.diff(ut, z)+ur*ut/r,
             ur*s.diff(uz, r)+uz*s.diff(uz, z)+s.diff(p, z))
    ledger.check("full cylindrical convective and pressure Euler residual vanishes",
                 all(reduced(value) == 0 for value in euler))
    ledger.check("axis regularity follows from finite even-profile limits",
                 s.limit(j1/r, r, 0) == kap/2 and s.limit(j0, r, 0) == 1)
    # At a simple J0 zero, substitute only that defining condition; the
    # Hessian and frequency are differentiated rather than assigned.
    hessian = s.hessian(psi, (r, z)).subs(z, 0)
    normalized = hessian/(amplitude*r*j1)
    actual = normalized.applyfunc(lambda v: s.simplify(s.expand_func(v)).subs(j0, 0))
    ledger.check("core streamfunction extremum has the derived nondegenerate Hessian",
                 s.simplify(actual-s.diag(-kap**2, -k**2)) == s.zeros(2, 2))
    transverse = s.Matrix([ur, uz]).jacobian((r, z)).subs(z, 0)
    transverse = transverse.applyfunc(lambda v: s.simplify(s.expand_func(v)).subs(j0, 0))
    ledger.check("actual transverse core frequency is obtained from the velocity jet",
                 s.simplify(transverse.det()-(amplitude*kap*k*j1)**2) == 0)
    ledger.check("finite-radius throughflow shear is retained",
                 s.simplify(s.expand_func(s.diff(ut/r, r)).subs(z, 0).subs(j0, 0)
                            +2*amplitude*lam*j1/r**2) == 0)

    a, b = s.symbols("a b", real=True)
    denominator = 1+a*a+b*b
    q = s.Matrix([2*a, 2*b, 1-a*a-b*b])/denominator
    ez = s.Matrix([0, 0, 1])
    coeff = ez-q*q[2]+s.I*q.cross(ez)
    ledger.check("rational stereographic directions remain on the exact unit sphere",
                 s.factor(q.dot(q)-1) == 0)
    ledger.check("helicity coefficient is exactly transverse and positive-curl",
                 s.simplify(q.dot(coeff)) == 0
                 and s.simplify(s.I*q.cross(coeff)-coeff) == s.zeros(3, 1))
    ledger.check("global coefficient bound follows from its actual squared norm",
                 s.simplify(s.conjugate(coeff).dot(coeff)-2*(1-q[2]**2)) == 0)

    # Three antipodal pairs with rational directions and positive weights
    # give one actual periodic field; every nonlinear Fourier mode is kept.
    wave_field = ({}, {}, {})
    for aa, bb in ((s.Rational(1, 3), s.Rational(1, 2)),
                   (s.Rational(1, 4), -s.Rational(2, 3)),
                   (s.Rational(3, 2), s.Rational(1, 5))):
        wave = tuple(q.subs({a: aa, b: bb}))
        vector = coeff.subs({a: aa, b: bb})/6
        for i in range(3):
            wave_field[i][wave] = vector[i]
            wave_field[i][tuple(-v for v in wave)] = s.conjugate(vector[i])
    actual_curl = fourier.curl(wave_field)
    ledger.check("commensurate paired field has exactly the same curl eigenvalue",
                 not fourier.divergence(wave_field)
                 and all(not fourier.add(actual_curl[i], fourier.scale(wave_field[i], -1))
                         for i in range(3)))
    pressure = fourier.scale(fourier.add(*(fourier.mul(v, v) for v in wave_field)),
                             -s.Rational(1, 2))
    convection = fourier.transport(wave_field, wave_field)
    ledger.check("the full periodic nonlinear Euler residual keeps all cross modes",
                 all(not fourier.add(convection[i], fourier.derivative(pressure, i))
                     for i in range(3)))
    wrong_curl = tuple(fourier.add(actual_curl[i], fourier.scale(wave_field[i], -2))
                       for i in range(3))
    ledger.check("a wrong curl normalization is exposed by the actual field",
                 any(wrong_curl))
    print("Scope: exact global toroidal field and periodic helicity identities;")
    print("tube persistence/approximation are analytic; local optical margins remain upstream")
    raise SystemExit(ledger.finish())


if __name__ == "__main__":
    main()
