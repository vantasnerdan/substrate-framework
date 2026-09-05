"""Exact Euler sideband algebra and the averaged-frequency counterexample."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0121-rotating-sideband")
    px, py, pz = s.symbols("px py pz", real=True)
    omega, carrier = s.symbols("Omega N", positive=True)
    p = s.Matrix([px, py, pz])
    norm = p.dot(p)
    projector = s.eye(3)-p*p.T/norm
    rotation_generator = omega*s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    generator = -2*projector*rotation_generator*projector
    frequency2 = 4*omega**2*pz**2/norm
    checks.check("complete pressure amplitude preserves the transverse plane",
                 s.simplify(p.T*generator) == s.zeros(1, 3))
    checks.check("co-rotating Euler amplitude has the derived squared frequency",
                 s.simplify(generator**2+frequency2*projector) == s.zeros(3))
    wrong = -projector*rotation_generator*projector
    checks.check("omitting frame rotation fails the actual amplitude polynomial",
                 s.simplify(wrong**2+frequency2*projector) != s.zeros(3))

    gradient = s.Matrix([s.diff(frequency2, x) for x in p])
    hessian = s.hessian(frequency2, tuple(p))
    c = s.Symbol("c", real=True)
    carrier_point = {px: carrier*s.sqrt(1-c**2), py: 0, pz: carrier*c}
    trace = s.simplify(s.trace(hessian).subs(carrier_point))
    gradient_norm = s.simplify(gradient.dot(gradient).subs(carrier_point))
    checks.check("sideband Hessian trace is derived before orientation averaging",
                 s.simplify(trace-8*omega**2*(1-3*c**2)/carrier**2) == 0)
    checks.check("orientation-induced first-jet variance is independently derived",
                 s.simplify(gradient_norm-64*omega**4*c**2*(1-c**2)/carrier**2) == 0)
    checks.check("oblique half-angle cosine gives positive scalar averaged curvature",
                 s.simplify((trace/6).subs(c, s.Rational(1, 2))-omega**2/(3*carrier**2)) == 0)
    checks.check("the same oblique candidate has a strictly nonzero omitted variance",
                 s.simplify((gradient_norm/3).subs(c, s.Rational(1, 2))-4*omega**4/carrier**2) == 0)
    checks.check("axial carrier removes the linear splitting but has negative curvature",
                 gradient_norm.subs(c, 1) == 0
                 and s.simplify((trace/6).subs(c, 1)+8*omega**2/(3*carrier**2)) == 0)

    f0, perturbation, linear, quadratic = s.symbols("f0 epsilon a b", real=True)
    plus = f0+perturbation*linear+perturbation**2*quadratic
    minus = f0-perturbation*linear+perturbation**2*quadratic
    exact_fourth = (plus**2+minus**2)/2
    copied_oscillator_fourth = ((plus+minus)/2)**2
    checks.check("actual initial fourth derivative exposes averaged-frequency closure",
                 s.expand(exact_fourth-copied_oscillator_fourth) == perturbation**2*linear**2)

    k = s.Matrix(s.symbols("k:3", real=True))
    # Apply the independently normalized fourth spherical moment to each
    # indexed contraction, rather than insert the tensor claimed in (5).
    identity = s.eye(3)
    fourth_contraction = s.Matrix(3, 3, lambda i, j: sum(
        k[a]*k[b]*(identity[a, b]*identity[i, j]
                   +identity[a, i]*identity[b, j]+identity[a, j]*identity[b, i])/15
        for a in range(3) for b in range(3)))
    screened = k.dot(k)*identity/3+fourth_contraction
    checks.check("vector screening contraction retains its longitudinal tensor",
                 s.simplify(s.Rational(3, 2)*screened
                            -(3*k.dot(k)*identity+k*k.T)/5) == s.zeros(3))
    print("Hessian trace:", trace)
    print("First-jet variance per unit slow |K|^2:", s.factor(gradient_norm/3))
    print("Scope: rotating-Euler comparison and closure counterexample, not EPS vector moduli")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
