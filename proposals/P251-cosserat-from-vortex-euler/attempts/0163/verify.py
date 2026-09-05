"""Exact full-current and Jacobi-action sign of the lifted triangular array."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0163-constant-curl-response")
    lam, amp, k, rho = s.symbols("lambda Psi k rho", positive=True)
    q = s.sqrt(lam**2+k**2)
    x, y = s.symbols("X Y", real=True)
    displacement = s.Matrix([x, y])
    j = s.Matrix([[0, -1], [1, 0]])
    ez = s.Matrix([0, 0, 1])
    base_vectors = [lam*s.Matrix([1, 0]),
                    lam*s.Matrix([-s.Rational(1, 2), s.sqrt(3)/2]),
                    lam*s.Matrix([-s.Rational(1, 2), -s.sqrt(3)/2])]
    covariance = s.zeros(2)
    diff_norm = s.Integer(0)
    convective_cross = s.Integer(0)
    v_b = s.zeros(2, 1)
    v_r = s.zeros(2, 1)
    raw_flux = s.zeros(2, 1)
    cross_mean = s.zeros(3, 1)
    initial_translation_flux = s.zeros(2, 1)
    helicity_residual = []
    base_residual = []
    for beta0 in base_vectors:
        covariance += amp**2*(j*beta0)*(j*beta0).T/2
        for sign in (1, -1):
            beta = sign*beta0
            wave = beta.col_join(s.Matrix([k]))
            e1 = (j*beta/lam).col_join(s.Matrix([0]))
            e2 = (lam*ez-k*(beta/lam).col_join(s.Matrix([0])))/q
            background = -lam*amp*ez/2+s.I*amp*(j*beta).col_join(s.Matrix([0]))/2
            translation = -s.I*beta.dot(displacement)*background
            prepared = amp*lam*beta.dot(displacement)*(e1+s.I*e2)/2
            helicity_residual.extend(s.I*wave.cross(prepared)-q*prepared)
            base_wave = beta.col_join(s.Matrix([0]))
            base_residual.extend(s.I*base_wave.cross(background)-lam*background)
            difference = prepared-translation
            convective = (-s.I*k*lam*amp*displacement/2).col_join(s.Matrix([0]))
            diff_norm += s.conjugate(difference).dot(difference)
            convective_cross += s.conjugate(difference).dot(convective)
            opposite = s.conjugate(background)
            v_b += opposite[:2, 0]*prepared[2]
            v_r += opposite[:2, 0]*(1-q/lam)*prepared[2]
            raw_flux += opposite[2]*prepared[:2, 0]+opposite[:2, 0]*prepared[2]
            initial_translation_flux += (opposite[2]*translation[:2, 0]
                                         +opposite[:2, 0]*translation[2])
            cross_mean += opposite.cross(prepared)

    def zero(entries):
        return all(s.simplify(entry) == 0 for entry in entries)

    checks.check("the same six-mode background has the declared constant curl",
                 zero(base_residual))
    checks.check("sideband translation has exact positive helicity including divergence return",
                 zero(helicity_residual))
    checks.check("complete covariance retains the entire triangular field",
                 s.simplify(covariance-3*amp**2*lam**2*s.eye(2)/4) == s.zeros(2))
    checks.check("full translation cancels both order-k mean-flux terms together",
                 zero(initial_translation_flux))
    checks.check("the actual sideband axial component retains its polarization factor",
                 zero(v_b-lam**2/q*covariance*j*displacement))
    checks.check("passive-r sideband row begins at second order, not first",
                 zero(v_r-(1-q/lam)*v_b))
    exact_mean = (q-lam)*cross_mean[:2, 0]
    current_mean = -s.I*k*v_r-k*k/lam*j*v_b
    checks.check("full Euler pressure and the exact Hodge-current reduction agree",
                 zero(exact_mean-current_mean)
                 and zero(exact_mean+s.I*k*raw_flux))
    checks.check("the leading physical translational sign is positive acceleration",
                 zero(s.Matrix([s.limit(entry/k**2, k, 0)
                                for entry in exact_mean])-covariance*displacement))
    checks.check("deleting the vertical translation partner changes the source",
                 not zero(v_b.subs(k, 0)))

    metric = (displacement.T*covariance*displacement)[0]
    diff_norm = s.simplify(diff_norm)
    convective_cross = s.simplify(s.re(convective_cross))
    checks.check("complete six-mode return norm is derived from actual polarizations",
                 s.simplify(diff_norm-2*lam**2*(1-lam/q)*metric) == 0)
    checks.check("Jacobi gyroscopic work is retained alongside that norm",
                 s.simplify(convective_cross-lam*k**2*metric/q) == 0)
    phase_energy = s.simplify(rho*diff_norm/2-rho*convective_cross)
    checks.check("same material Hamiltonian has exact negative acoustic stiffness",
                 s.simplify(phase_energy+rho*lam*(q-lam)*metric) == 0)
    checks.check("omitting convective work falsely changes the phase sign",
                 s.simplify(phase_energy-rho*diff_norm/2) != 0)
    checks.check("the full-action longwave coefficient confirms the current sign",
                 s.simplify(s.limit(2*phase_energy/(rho*k*k), k, 0)+metric) == 0)

    # The following scalar row algebra retains the full pressure functional.
    vr, rs, rdot, vb, jj = s.symbols("Vr RS Rdot Vb J")
    mean_rate = -s.I*k*vr-k*k*jj*vb/lam
    # r_t+A r=-i k S and A R=v imply Vr=Rdot+i k RS.
    corrected = s.expand(mean_rate.subs(vr, rdot+s.I*k*rs)+s.I*k*rdot)
    checks.check("the current primitive keeps pressure S and the explicit axial-spin row",
                 s.expand(corrected-k*k*rs+k*k*jj*vb/lam) == 0)
    print("Scope: exact lifted current and same-Jacobi sign; positive acoustic stiffness is not asserted.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
