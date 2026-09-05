"""Exact rational repair of the primary Table 2 displacement, not MHD dynamics."""

import importlib.util
from pathlib import Path

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    source = Path(__file__).resolve().parents[1]/"0040"/"fourier_orbit.py"
    spec = importlib.util.spec_from_file_location("p251_exact_fourier", source)
    fourier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fourier)
    ledger = CheckLedger("P251-0089-rational-material-witness")
    # Integers encode the four-decimal Table 2 coefficients exactly.
    rows = [
        ((1, 1, 0), (-1007, 1007, 5), (-1007, 1007, -6)),
        ((1, -1, 0), (-600, -600, -3), (-600, -600, 4)),
        ((1, 3, 0), (-108, 36, 1), (-108, 36, 0)),
        ((1, -3, 0), (-181, -60, -1), (-181, -60, 1)),
        ((3, 1, 0), (-36, 108, -1), (36, -108, 0)),
        ((3, -1, 0), (-60, -181, 1), (60, 181, 1)),
        ((3, 3, 0), (-139, 139, 0), (139, -139, 0)),
        ((3, -3, 0), (-83, -83, 0), (83, 83, 0)),
    ]
    raw = ({}, {}, {})
    for k, real, imag in rows:
        for j in range(3):
            value = (s.Integer(real[j])+s.I*imag[j])/10000
            raw[j][k] = value
            raw[j][tuple(-a for a in k)] = s.conjugate(value)
    xi = fourier.leray(raw)
    ledger.check("printed rounding is explicitly detected before repairing incompressibility",
                 fourier.divergence(raw) != {})
    ledger.check("rational projected displacement is exactly divergence free",
                 fourier.divergence(xi) == {})
    ledger.check("real smooth finite Fourier displacement retains conjugate coefficients",
                 all(s.conjugate(value) == component[tuple(-a for a in k)]
                     for component in xi for k, value in component.items()))
    u = (fourier.scale(fourier.trig(1, 2, "sin"), -1),
         fourier.trig(0, 2),
         fourier.add(fourier.trig(1, 2), fourier.scale(fourier.trig(0, 2, "sin"), -1)))
    ledger.check("source-convention two-wave field has positive Beltrami eigenvalue two",
                 fourier.curl(u) == tuple(fourier.scale(a, 2) for a in u)
                 and fourier.divergence(u) == {})
    f = fourier.cross(xi, u)
    w = fourier.curl(f)
    helicity = fourier.inner(f, w)
    curl_norm = fourier.inner(w, w)
    stiffness = s.factor(2*helicity-curl_norm)
    mass = fourier.inner(xi, xi)
    print("EXACT mean K/rho:", stiffness)
    print("EXACT mean material mass/rho:", mass)
    print("EXACT mean F.curl F:", helicity)
    print("EXACT mean |curl F|²:", curl_norm)
    print("EXACT K/mass:", s.factor(stiffness/mass))
    ledger.check("actual rational trial has strictly positive material Jacobi stiffness",
                 stiffness > 0 and mass > 0)
    # Independent form: pressure Hessian minus full background convective norm.
    pressure = fourier.scale(fourier.add(*(fourier.mul(a, a) for a in u)), -s.Rational(1, 2))
    transport = tuple(fourier.add(*(fourier.mul(u[j], fourier.derivative(xi[i], j))
                                    for j in range(3))) for i in range(3))
    pressure_hessian = sum(fourier.mul(fourier.mul(xi[i], xi[j]),
                           fourier.derivative(fourier.derivative(pressure, i), j))
                           .get((0, 0, 0), 0) for i in range(3) for j in range(3))
    direct = s.expand(pressure_hessian-fourier.inner(transport, transport))
    print("EXACT direct pressure-Hessian minus convective norm:", direct)
    ledger.check("independent full Jacobi pressure-and-transport form gives the same rational sign",
                 direct == stiffness)
    ledger.check("omitting helicity or reversing its sign removes the positive witness",
                 -curl_norm < 0 and -2*helicity-curl_norm < 0)
    ledger.check("the imported MHD potential has exactly the opposite functional sign",
                 s.factor((curl_norm-2*helicity)/2) == -stiffness/2)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
