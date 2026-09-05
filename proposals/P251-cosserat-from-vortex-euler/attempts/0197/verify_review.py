"""Independent Cartesian phase/Haar audit of the exposed first-shell constants."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0197-independent-D-energy")
    kap = s.Matrix(s.symbols("kx ky kz", real=True))
    disp = s.Matrix(s.symbols("Dx Dy Dz", real=True))
    aa, bb, parameter = s.symbols("A B t", real=True)
    units = [s.eye(3)[:, i] for i in range(3)]
    d = kap[1]*disp[1]-kap[2]*disp[2]

    def average(poly):
        # Direct fourth moment of an orthonormal pair; no source helper.
        result = 0
        for powers, coeff in s.Poly(s.expand(poly), *kap, *disp).terms():
            ki = [i for i in range(3) for _ in range(powers[i])]
            di = [i for i in range(3) for _ in range(powers[i+3])]
            if len(ki) != 2 or len(di) != 2:
                raise ValueError("expected the actual quadratic strain polynomial")
            i, j = ki
            a, b = di
            moment = s.Rational(2, 15)*int(i == j and a == b)
            moment -= s.Rational(1, 30)*(int(i == a and j == b)+int(i == b and j == a))
            result += coeff*moment
        return s.factor(result)

    norm_after = 0
    current_after = 0
    base_ad_norm = 0
    projected_current_norm = 0
    omitted_axial_energy = 0
    for axis, amp, sin_unit, axial_sign in (
        (1, bb, -units[2], -1), (2, aa, units[1], 1)
    ):
        normal = units[axis]
        uc, us = amp*units[0], amp*sin_unit
        checks.check(
            f"wave {axis} has the actual negative helicity in Cartesian phases",
            normal.cross(us) == -uc and -normal.cross(uc) == -us,
        )
        qc = disp[axis]*kap.cross(us)
        qs = -disp[axis]*kap.cross(uc)
        bc0, bs = qc-kap.dot(uc)*disp, qs-kap.dot(us)*disp
        bc = bc0+axial_sign*d*amp*units[0]/2
        cc = kap.dot(uc)*disp+disp.dot(uc)*kap
        cs = kap.dot(us)*disp+disp.dot(us)*kap
        p = s.eye(3)-normal*normal.T

        def minus(vc, vs):
            return (p*vc-normal.cross(vs))/2, (p*vs+normal.cross(vc))/2

        mc, ms = minus(bc, bs)
        nc, ns = minus(cc, cs)
        fc, fs = bc-mc+parameter*nc, bs-ms+parameter*ns
        checks.check(
            f"wave {axis} return is a genuine solenoidal stationary helicity field",
            s.simplify(normal.dot(-mc+parameter*nc)) == 0
            and s.simplify(normal.cross(-ms+parameter*ns)+(-mc+parameter*nc)) == s.zeros(3, 1),
        )
        norm_after += (fc.dot(fc)+fs.dot(fs))/2
        current_after += (cc.dot(fc)+cs.dot(fs))/2
        base_ad_norm += (kap.dot(uc)**2+kap.dot(us)**2)*disp.dot(disp)/2
        projected_current_norm += (nc.dot(nc)+ns.dot(ns))/2
        omitted_c, omitted_s = minus(bc0, bs)
        bad_c, bad_s = bc0-omitted_c+parameter*nc, bs-omitted_s+parameter*ns
        omitted_axial_energy += (bad_c.dot(bad_c)+bad_s.dot(bad_s))/2

    total_energy = aa**2+bb**2
    actual_h = average(norm_after-base_ad_norm)
    actual_r = average(base_ad_norm+current_after)
    checks.check(
        "Cartesian full-shell energy independently gives the 47 over 240 constant",
        s.factor(actual_h-total_energy*(parameter**2/s.Integer(15)-s.Rational(47, 240))) == 0,
    )
    checks.check(
        "independent physical current gives the 13 over 120 constant",
        s.factor(actual_r-total_energy*(8*parameter+13)/120) == 0,
    )
    checks.check(
        "complete correlated current projector has the stated positive norm",
        s.factor(average(projected_current_norm)-total_energy/15) == 0,
    )
    checks.check(
        "omitting the actual stationary axial return changes the exposed energy",
        s.factor(average(omitted_axial_energy-base_ad_norm)-actual_h) != 0,
    )
    range_energy, restoring = s.symbols("R a", real=True)
    selected_t = s.solve(actual_r+restoring, parameter)[0]
    mismatch = s.factor((actual_h+range_energy+actual_r).subs(parameter, selected_t))
    expected = range_energy+(120*restoring-total_energy)*(120*restoring+19*total_energy)/(960*total_energy)
    checks.check(
        "actual energy minus independently observed stiffness derives the stated quadratic",
        s.factor(mismatch-expected) == 0,
    )
    print("Independent Cartesian first-shell energy:", actual_h)
    print("Independent physical mean acceleration:", actual_r)
    print("This check does not supply the unresolved V-column or a parent continuum.")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
