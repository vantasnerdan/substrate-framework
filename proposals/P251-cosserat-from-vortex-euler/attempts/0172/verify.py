"""Exact exposing checks of the prepared-history action and physical normal form."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_phase import physical_configuration_chart
from substrate_framework.micropolar import MicropolarCoefficients, micropolar_fourier_stiffness
from substrate_framework.verification import CheckLedger


def main() -> int:
    checks = CheckLedger("P251-0172")
    t, w, a, b, c = sp.symbols("t omega A B c", positive=True)
    lam = sp.symbols("lambda", positive=True)
    f = sp.cos(sp.sqrt(lam) * t)
    g = sp.sin(sp.sqrt(lam) * t) / sp.sqrt(lam)
    f2 = (b * f.diff(lam) + a**2 * f.diff(lam, 2) / 2).subs(lam, w**2)
    g2 = (b * g.diff(lam) + a**2 * g.diff(lam, 2) / 2).subs(lam, w**2)
    f, g = f.subs(lam, w**2), g.subs(lam, w**2)
    def simplify(value):
        return sp.trigsimp(sp.simplify(sp.expand_trig(value)))
    w2 = simplify(f2 * g.diff(t) + f * g2.diff(t) - g2 * f.diff(t) - g * f2.diff(t))
    n2 = simplify(f2.diff(t) * g.diff(t, 2) + f.diff(t) * g2.diff(t, 2)
                  - g2.diff(t) * f.diff(t, 2) - g.diff(t) * f2.diff(t, 2))
    hfun = (w**2 * t**2 - sp.sin(w * t)**2) / (4 * w**4)
    qfun = sp.sin(w * t)**2 / (2 * w**2)
    checks.check("derived coherent Wronskian retains the group variance", simplify(w2 + a**2 * hfun) == 0)
    checks.check("derived coherent stiffness retains both moving terms", simplify(n2 - 2 * w**2 * w2 - b - a**2 * (w**2 * hfun - qfun)) == 0)
    checks.check("ratio curvature is not the average squared-frequency curvature", simplify(n2 - w**2 * w2 - b) != 0)
    checks.check("actual mass-rate connection agrees with retained memory", simplify(hfun.diff(t) + (sp.sin(w*t)*sp.cos(w*t)-w*t)/(2*w**3)) == 0)
    checks.check("group-zero removes this variance, not other physical connections", w2.subs(a, 0) == 0 and simplify((n2-w**2*w2).subs(a, 0)-b) == 0)
    r = sp.symbols("r", real=True)
    for name, free in (("angle", f), ("rate", g)):
        memory = sp.integrate(sp.sin(w*(t-r))/w * free.subs(t, r), (r, 0, t))
        row = qfun * free + (sp.sin(w*t)*sp.cos(w*t)-w*t)/(2*w**3) * free.diff(t)
        checks.check(f"actual retained convolution {name} column", simplify(memory-row) == 0)

    z = sp.symbols("z", real=True)
    transverse = 3 * sp.integrate((1-z**2)*z**2, (z,-1,1))/4
    longitudinal = 3 * sp.integrate(z**4, (z,-1,1))/2
    checks.check("one common laboratory K gives fourth-Haar transverse factor", transverse == sp.Rational(1,5))
    checks.check("one common laboratory K gives fourth-Haar longitudinal factor", longitudinal == sp.Rational(3,5))
    checks.check("odd common-K sideband tensor averages to zero", sp.integrate(z**3, (z,-1,1)) == 0)
    delta, omega, p, ap, nn = sp.symbols("delta Omega p a_packet N", positive=True)
    cd = sp.sqrt(2)*omega*delta
    slope = -2*ap*omega*cd/(nn*p)
    curvature = 3*ap**2*omega*cd/(2*nn*p**2)
    checks.check("packet-own group error has the small delta scale", sp.simplify(slope**2/(4*omega**2*curvature)-2*sp.sqrt(2)*delta/(3*nn)) == 0)

    k, rho, j, inertia, ma, j2, speed, bt, gap = sp.symbols("k rho j I ma j2 a BT gap", positive=True)
    def jet(matrix):
        return matrix.applyfunc(lambda entry: sp.series(entry, k, 0, 3).removeO().expand())
    physical = sp.Matrix([[1-inertia*k**2/(4*rho), -j*k/(2*rho)], [k/2, 1]])
    inverse = jet(physical.inv())
    raw_mass = sp.diag(rho+ma*k**2, j+j2*k**2)
    raw_stiff = sp.diag(rho*speed*k**2, j*gap+(j*bt+gap*j2)*k**2)
    mass = jet(inverse.T*raw_mass*inverse)
    stiff = jet(inverse.T*raw_stiff*inverse)
    expected_mass = sp.diag(rho+(ma+inertia/2-j/4)*k**2, j+(j2-j**2/(4*rho))*k**2)
    expected_stiff = sp.Matrix([[(rho*speed+j*gap/4)*k**2,-j*gap*k/2],[-j*gap*k/2,j*gap+(j*bt+gap*j2-j**2*gap/(2*rho))*k**2]])
    checks.check("literal acoustic I and optical j yield the complete physical mass", sp.simplify(mass-expected_mass) == sp.zeros(2))
    checks.check("the same physical map pulls back the complete potential", sp.simplify(stiff-expected_stiff) == sp.zeros(2))
    normalizer = sp.eye(2)-sp.diag(1/rho,1/j)*(mass-sp.diag(rho,j))/2
    norm_mass, norm_stiff = jet(normalizer.T*mass*normalizer), jet(normalizer.T*stiff*normalizer)
    ct = j*bt-j**2*gap/(4*rho)
    checks.check("derived derivative map normalizes BOTH inertias", sp.simplify(norm_mass-sp.diag(rho,j)) == sp.zeros(2))
    checks.check("derived transverse curvature includes the small-inertia subtraction", sp.simplify(norm_stiff[1,1]-j*gap-ct*k**2) == 0)
    checks.check("gradient tag inertia and carrier mass cancel from normalized curvature", not norm_stiff.has(inertia, ma, j2))
    checks.check("positive optical frequency curvature alone is not a sign proof", ct.subs(j, 8*rho*bt/gap) < 0)
    checks.check("finite positive density can make the coupled margin strict", sp.simplify(ct.subs(j, rho*bt/gap)-3*rho*bt**2/(4*gap)) == 0)
    bl = sp.symbols("BL", positive=True)
    cl = j*bl
    coefficients = MicropolarCoefficients(0, rho*speed, j*gap/4, (cl-ct)/2, ct/2, ct/2)
    canonical = sp.Matrix(micropolar_fourier_stiffness([0,0,k], coefficients))
    helicity = sp.Matrix([1, sp.I, 0])/sp.sqrt(2)
    embedding = sp.zeros(6,2)
    embedding[:3,0], embedding[3:,1] = helicity, helicity
    checks.check("canonical micropolar API reproduces the entire transverse normal pencil", sp.simplify(embedding.conjugate().T*canonical*embedding-norm_stiff) == sp.zeros(2))
    checks.check("canonical longitudinal spin curvature is separately retained", sp.simplify(canonical[5,5]-j*gap-cl*k**2) == 0)
    value = sp.symbols("value")
    determinant = (norm_stiff-value*sp.diag(rho,j)).det().expand()
    acoustic = sp.expand(determinant.subs(value, speed*k**2)).coeff(k,2)
    optical = sp.expand(determinant.subs(value, gap+bt*k**2)).coeff(k,2)
    checks.check("both physical dispersion branches survive the field normalization", acoustic == 0 and sp.simplify(optical) == 0)

    # The canonical configuration API supplies all moving terms, rather than
    # an independently copied phase elimination. A variable observer exposes them.
    phase = sp.Matrix([[0,1],[-1,0]])
    observer = sp.Matrix([[1+t**2, t]])
    chart = physical_configuration_chart(phase, sp.zeros(2), observer,
        configuration_rate=observer.diff(t), configuration_acceleration=observer.diff(t,2),
        generator_rate=sp.zeros(2), momentum=sp.Matrix([[0,1]]))
    checks.check("actual chart mass is the inverse physical Wronskian", sp.simplify(chart.mass[0,0]-1/(1-t**2)) == 0)
    checks.check("actual moving chart retains a measured/canonical momentum difference", chart.momentum_difference != sp.zeros(1,2))
    clock, amplitude = sp.Function("gamma")(t), sp.Function("c")(t)
    moving_mass = -1/(clock*amplitude**2)
    ell = amplitude.diff(t)/amplitude
    connection = moving_mass*(clock.diff(t)/clock+ell)
    checks.check("integrated physical spin retains the exact clock-current connection", sp.simplify(-moving_mass*ell-moving_mass.diff(t)-connection) == 0)
    # A scalar jet probe exposes both the actual initial dipole and its current.
    angle = 1+t+t**2
    target = sp.diff(moving_mass*angle, t)+connection*angle
    checks.check("differentiated hybrid current equals the actual matched spin row", sp.simplify(target-moving_mass*(angle.diff(t)-ell*angle)) == 0)
    x, y, zz = sp.symbols("x y z", real=True)
    coords = (x,y,zz)
    vector = sp.Matrix([x*y+zz**2, y*zz+x**2, zz*x+y**2])
    grad = vector.jacobian(coords)
    flux = vector*sp.trace(grad)-grad*vector
    checks.check("curvature representative change retains its exact surface flux", sp.expand(sp.trace(grad)**2-sp.trace(grad*grad)-sum(sp.diff(flux[i],coords[i]) for i in range(3))) == 0)
    checks.check("positive curvature representative has the required trace margin", sp.simplify((3*coefficients.trace_curvature+coefficients.symmetric_curvature)-(3*cl-2*ct)/2) == 0)
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
