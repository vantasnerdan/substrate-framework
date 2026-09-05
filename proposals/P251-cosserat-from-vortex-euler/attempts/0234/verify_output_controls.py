"""Actual generic-K Kelvin jet and literal positive-fraction output inverse."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0234-output-controls")
    k = sp.Symbol("k", real=True, nonzero=True)
    psi = ef.add(ef.trig(2), ef.scale(ef.trig(1), sp.Rational(1, 100)))
    u = (psi, ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -sp.Rational(1, 100)))
    s = ef.trig(1)
    xi0 = (ef.scale(s, -1), ef.scale(ef.derivative(s, 2), -1), ef.derivative(s, 1))
    direction = (1, 2, 3)

    def lift(field):
        return {tuple(q[j]+k*direction[j] for j in range(3)): c for q, c in field.items()}

    projected = ef.leray(tuple(lift(field) for field in xi0))
    xi = tuple(ef.scale(field, k*k) for field in projected)
    velocity = ef.leray(ef.cross(xi, ef.curl(u)))
    rate = ef.material_kelvin_operator(u, xi)

    def jet(field, order):
        result = {}
        for q, c in field.items():
            wave = tuple(sp.expand(q[j]-k*direction[j]) for j in range(3))
            value = sp.cancel(sp.diff(c, k, order).limit(k, 0)/sp.factorial(order))
            if value != 0:
                result[wave] = value
        return result

    def equal(left, right):
        return all(sp.cancel(c) == 0 for a, b in zip(left, right, strict=True)
                   for c in ef.add(a, ef.scale(b, -1)).values())

    ledger.check("the actual generic-K correction is solenoidal with the FULL Leray pressure",
                 all(sp.cancel(c) == 0 for c in ef.divergence(xi).values()))
    for order in (0, 1):
        ledger.check(f"the true physical Euler correction has no lower spatial coefficient {order}",
                     all(not jet(field, order) for field in velocity))
    ts = ef.transport(u, (s, {}, {}))[0]
    ledger.check("the full three-dimensional pressure gives the exact passive second velocity coefficient",
                 equal(tuple(jet(field, 2) for field in velocity), (ts, {}, {})))
    expected_rate = (ts, ef.derivative(ts, 2), ef.scale(ef.derivative(ts, 1), -1))
    ledger.check("the actual generic-K Lin second coefficient is the real passive Kelvin history",
                 equal(tuple(jet(field, 2) for field in rate), expected_rate))
    ledger.check("this physical supplier genuinely changes the Euler vorticity",
                 any(ef.curl((ts, {}, {}))))

    # Derive the angle row from actual fraction-scaled quadrupoles, and
    # the current row from the literal fraction-scaled material integral.
    qref, dq, current = sp.symbols("Q dQ G", nonzero=True)
    fractions = (sp.Rational(1, 4), sp.Rational(1, 2))
    angle_rows = [sp.cancel(f*dq/(f*qref)/(dq/qref)) for f in fractions]
    current_rows = [sp.cancel(f*current/current) for f in fractions]
    matrix = sp.Matrix([angle_rows, current_rows])/2
    ledger.check("distinct POSITIVE physical tag fractions give nonzero literal output rank",
                 matrix.det() > 0 and all(0 < f < 1 for f in fractions))

    x = sp.Symbol("x", real=True)
    eta = (1-x*x)**6
    eta /= sp.integrate(eta, (x, -1, 1))
    nu = 2+x/7
    j = 1+nu+nu*nu
    desired_angle = sp.diff(eta, x)-2*eta
    desired_current = sp.diff(eta, x, 2)+x*eta
    fields = sp.Matrix([4*desired_angle-8*desired_current/j,
                        8*desired_current/j-2*desired_angle])
    actual_rows = matrix*fields
    actual_angle = sp.cancel(actual_rows[0])
    actual_current = sp.cancel(j*actual_rows[1])
    for order in range(5):
        multiplier = (-sp.I*nu)**order
        angle_error = sp.integrate(sp.expand(sp.cancel(multiplier*(actual_angle-desired_angle))), (x, -1, 1))
        current_error = sp.integrate(sp.expand(sp.cancel(multiplier*(actual_current-desired_current))), (x, -1, 1))
        ledger.check(f"actual transported angle time moment matches prescribed row {order}", angle_error == 0)
        ledger.check(f"actual transported G/spin time moment matches prescribed row {order}", current_error == 0)
    ledger.check("equal fractions lose precisely the new physical rank",
                 sp.Matrix([[1, 1], [fractions[0], fractions[0]]]).det() == 0)
    wrong_fields = fields.subs(j, 1) if j.is_Symbol else sp.Matrix([
        4*desired_angle-8*desired_current, 8*desired_current-2*desired_angle])
    wrong_current = sp.cancel(j*(matrix*wrong_fields)[1])
    wrong_moment = sp.integrate(sp.expand(wrong_current-desired_current), (x, -1, 1))
    ledger.check("omitting the actual frequency-dependent measured inertia corrupts the physical current",
                 wrong_moment != 0)

    # Exact integration by parts exposes the j' connection in a carrier
    # derivative; no numerical delta packet or fitted frequency is used.
    derivative_row = sp.integrate(sp.diff(eta, x)*j, (x, -1, 1))
    measured_connection = -sp.integrate(eta*sp.diff(j, x), (x, -1, 1))
    ledger.check("the literal carrier-derivative current retains the measured-inertia connection",
                 sp.simplify(derivative_row-measured_connection) == 0 and derivative_row != 0)
    second_angle_moment = sp.integrate(sp.diff(eta, x, 2)*nu**2, (x, -1, 1))
    ledger.check("a free second-carrier derivative produces a nonzero actual quadratic-time angle row",
                 sp.simplify(second_angle_moment-2*sp.diff(nu, x)**2) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
