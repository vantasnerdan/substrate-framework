"""Exact local algebra supporting the supervisor's 0037 mathematical notes.

The continuum norm estimates and cutoff/domain arguments are in the notes.
These identities expose coordinate, curl, scaling and derivative-placement
errors; they do not verify the unfinished Euler scattering theorem.
"""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P253-0037-supervisor-analysis")
    r, eps, s, W = sp.symbols("r epsilon s W", positive=True)
    z, t = sp.symbols("z t", real=True)
    F, G = sp.Function("F"), sp.Function("G")
    radial = eps**2 * F(r / eps)
    g = G(z).diff(z)
    psi = radial * g
    vr, vz = -psi.diff(z) / r, psi.diff(r) / r
    omega = vr.diff(z) - vz.diff(r)
    eta = omega / r
    boundary_density = -(psi.diff(r) / r).diff(r) - psi.diff(z, 2) / r
    ledger.check("S2 exact radial flux and axis-boundary sign",
                 sp.simplify(r * eta - boundary_density) == 0)
    ledger.check("S2 reversed axial-Laplacian sign is exposed",
                 sp.simplify(r * eta + (psi.diff(r) / r).diff(r)
                             - psi.diff(z, 2) / r) != 0)

    primitive_displacement = -(1 / (W * r)) * (
        (radial.diff(r, 2) - radial.diff(r) / r) * G(z)
        + radial * G(z).diff(z, 2))
    ledger.check("S3-S4 exact compact-generator curl identity",
                 sp.simplify(W * primitive_displacement.diff(z) - omega) == 0)

    changed_measure_density = ((vr**2 + vz**2) * r * eps).subs(r, eps * s).doit()
    expected_density = (eps**2 * F(s).diff(s)**2 * g**2 / s
                        + eps**4 * F(s)**2 * g.diff(z)**2 / s)
    ledger.check("S1 exact radial energy scaling including measure",
                 sp.simplify(changed_measure_density - expected_density) == 0)

    a, b, length = sp.symbols("a b length", positive=True)
    B, P = sp.Function("B"), sp.Function("P")
    d = a - b  # The note assumes d != 0.
    cutoff = B((z - a * t) / length)
    profile = P(z - d * t)
    integrand = cutoff.diff(z) * profile + cutoff * profile.diff(z)
    source_derivative = (cutoff.diff(z) * profile.diff(z)
                         + cutoff * profile.diff(z, 2))
    ledger.check("A5 exact transport homological identity",
                 sp.simplify(a * integrand.diff(z) - b * source_derivative
                             + integrand.diff(t)) == 0)
    reduced = -(cutoff * profile).diff(t) / d - b * cutoff.diff(z) * profile / d
    ledger.check("A6 output derivative yields the source-speed factor",
                 sp.simplify(integrand - reduced) == 0)
    wrong_speed = -(cutoff * profile).diff(t) / d - a * cutoff.diff(z) * profile / d
    ledger.check("A6 exchanging source and target speeds is exposed",
                 sp.simplify(integrand - wrong_speed) != 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
