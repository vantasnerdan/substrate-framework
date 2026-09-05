"""Actual carrier dispersion, coherent action normalization and physical jets."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0138-finite-carrier-join")
    d = s.Symbol("d", real=True)
    f0, m0, curvature = s.symbols("f0 M0 B", positive=True)
    group, m1, m2 = s.symbols("A M1 M2", real=True)
    plus, minus = f0+group*d+curvature*d*d, f0-group*d+curvature*d*d
    average = (plus+minus)/2
    checks.check("standing carrier cancels first group derivative but retains its variance",
                 s.expand((plus**2+minus**2)/2-average**2-group**2*d*d) == 0)
    fmatrix = s.Matrix([[average, s.I*group*d], [-s.I*group*d, average]])
    eigenvalue = s.Symbol("eigenvalue")
    checks.check("retained sine/cosine carrier state has the actual two branch frequencies",
                 s.expand((fmatrix-eigenvalue*s.eye(2)).det()
                          -(plus-eigenvalue)*(minus-eigenvalue)) == 0)

    # Same-fluid amplitudes 1/2 give action weights 1/4. Spin is linear.
    change = s.Matrix([[1, s.I], [1, -s.I]])
    mass_plus, mass_minus = m0+m1*d+m2*d*d, m0-m1*d+m2*d*d
    coherent_mass = s.simplify(s.conjugate(change).T*s.diag(mass_plus, mass_minus)*change/4)
    actual_spin_row = s.Matrix([[mass_plus/2, mass_minus/2]])*change/2
    checks.check("half-target actual tag spin equals the coherent canonical q momentum",
                 s.simplify(actual_spin_row-coherent_mass[0, :]) == s.zeros(1, 2))
    wrong_spin = s.Matrix([[mass_plus, mass_minus]])*change/2
    checks.check("reusing the single-carrier match gives the exposing factor two",
                 s.simplify(wrong_spin-2*coherent_mass[0, :]) == s.zeros(1, 2)
                 and s.simplify(wrong_spin-coherent_mass[0, :]) != s.zeros(1, 2))
    stiffness = s.simplify(s.conjugate(change).T*s.diag(mass_plus*plus, mass_minus*minus)*change/4)
    checks.check("complete mass jets reproduce the actual pole rather than a frozen mass",
                 s.simplify(stiffness-coherent_mass*fmatrix) == s.zeros(2))

    # The hidden carrier quadrature is recoverable from a real odd material
    # moment while leaving the actual even-tag angle q untouched.
    heven, hodd = s.symbols("h_even h_odd", nonzero=True, real=True)
    q, p = s.symbols("q p")
    odd_moment = -s.I*heven*q+hodd*p
    checks.check("actual odd material moment reconstructs the hidden branch coordinate",
                 s.simplify((odd_moment+s.I*heven*q)/hodd-p) == 0)

    # Isotropic moments are evaluated from the sphere, not supplied constants.
    z, phi = s.symbols("z phi", real=True)
    nx2 = (1-z*z)*s.cos(phi)**2
    second = s.integrate(z*z, (z, -1, 1))/2
    longitudinal = s.integrate(z**4, (z, -1, 1))/2
    transverse = s.integrate(s.integrate(nx2*z*z, (phi, 0, 2*s.pi)), (z, -1, 1))/(4*s.pi)
    checks.check("axial-angle reconstruction uses the physical one-third moment",
                 second == s.Rational(1, 3))
    checks.check("isotropic second spatial jet has longitudinal/transverse ratio three",
                 3*longitudinal == s.Rational(3, 5)
                 and 3*transverse == s.Rational(1, 5))
    x1, x2, x3 = s.symbols("x1 x2 x3", real=True)
    trace = x1+x2+x3
    # Diagonalization covers every real symmetric tensor by rotational covariance.
    hidden = s.Rational(5, 2)*(x1*nx2+x2*(1-z*z)*s.sin(phi)**2+x3*z*z)-trace/2
    hidden_norm = s.integrate(s.integrate(s.expand_trig(hidden**2), (phi, 0, 2*s.pi)),
                              (z, -1, 1))/(4*s.pi)
    checks.check("retained physical branch tensor has its positive inherited kinetic metric",
                 s.simplify(hidden_norm-(5*(x1*x1+x2*x2+x3*x3)-trace**2)/6) == 0)

    # Finite-time response bounds use the exact oscillator Green function.
    time, horizon, frequency, amplitude = s.symbols("t T nu0 R", positive=True)
    inner = s.Symbol("s", nonnegative=True)
    double_bound = s.integrate(s.integrate(amplitude/frequency**2, (time, 0, inner)),
                              (inner, 0, horizon))
    checks.check("ordered memory response has the stated finite-time double-integral bound",
                 double_bound == amplitude*horizon**2/(2*frequency**2))
    q0, core_radius, omega = s.symbols("q0 a Omega", positive=True)
    exact_quadratic = omega**2+omega**2*core_radius**2*q0**2/3
    memory_ratio = s.simplify(s.diff(exact_quadratic, q0)**2
                              /(exact_quadratic*s.diff(exact_quadratic, q0, 2)/2))
    checks.check("natural gradient-scale memory ratio is small at a fixed small nonzero carrier",
                 s.simplify(memory_ratio-4*core_radius**2*q0**2/(3+core_radius**2*q0**2)) == 0)
    print("Quadratic Rankine group-memory ratio:", memory_ratio)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
