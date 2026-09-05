"""Full Euler/Lin material-label feedback and exact physical-current identities."""

import sympy as s

from substrate_framework import euler_fourier as f
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0150-material-label-feedback")
    k = s.Symbol("k", positive=True)
    psi = f.mul(f.trig(0), f.trig(1))
    base = (f.scale(f.derivative(psi, 1), -1), f.derivative(psi, 0), {})
    zeta = f.curl(base)[2]
    carrier = {(0, 0, k): s.Integer(1)}
    tag = f.mul(carrier, f.trig(0, kind="sin"))
    seed = (f.mul(carrier, f.trig(1)),
            f.mul(carrier, f.trig(0, kind="sin")),
            f.mul(carrier, f.mul(f.trig(0), f.trig(1, kind="sin"))))
    velocity = f.leray(seed)
    first = f.transport(base, velocity)
    second = f.transport(velocity, base)
    rate = f.leray(tuple(f.scale(f.add(first[j], second[j]), -1) for j in range(3)))
    scalar_transport_tag = f.transport(base, ({}, {}, tag))[2]
    tag_rate = f.add(velocity[2], f.scale(scalar_transport_tag, -1))
    vertical_vorticity = f.curl(velocity)[2]
    q = f.add(vertical_vorticity, f.scale(f.derivative(f.mul(zeta, tag), 2), -1))
    qrate = f.add(f.curl(rate)[2],
                 f.scale(f.derivative(f.mul(zeta, tag_rate), 2), -1))
    aq = f.transport(base, ({}, {}, q))[2]
    advected_zeta = f.transport(velocity, ({}, {}, zeta))[2]
    exact_label_residual = f.add(qrate, aq, advected_zeta)
    checks.check("full nonlocal Euler velocity satisfies the exact rearrangement-label equation",
                 all(s.simplify(value) == 0 for value in exact_label_residual.values()))
    wrong_qrate = f.add(f.curl(rate)[2], f.derivative(f.mul(zeta, tag_rate), 2))
    wrong_q = f.add(vertical_vorticity, f.derivative(f.mul(zeta, tag), 2))
    wrong_residual = f.add(wrong_qrate,
                           f.transport(base, ({}, {}, wrong_q))[2], advected_zeta)
    checks.check("a reversed material stretching sign is detected",
                 any(s.simplify(value) != 0 for value in wrong_residual.values()))

    def green(scalar):
        return {wave: value/(wave[0]**2+wave[1]**2)
                for wave, value in scalar.items() if wave[0] != 0 or wave[1] != 0}

    def inverse_curl(scalar):
        potential = green(scalar)
        return (f.derivative(potential, 1), f.scale(f.derivative(potential, 0), -1), {})

    def gradient_green(scalar):
        potential = green(scalar)
        return (f.derivative(potential, 0), f.derivative(potential, 1), {})

    kq = inverse_curl(q)
    stretching = inverse_curl(f.derivative(f.mul(zeta, tag), 2))
    divergence_return = gradient_green(f.derivative(velocity[2], 2))
    mean_wave = (0, 0, k)
    reconstructed = tuple(f.add(kq[j], stretching[j], divergence_return[j],
                                {mean_wave: velocity[j].get(mean_wave, 0)}) for j in range(2))
    checks.check("physical velocity retains both stretching and divergence returns",
                 all(s.simplify(value) == 0
                     for j in range(2)
                     for value in f.add(reconstructed[j], f.scale(velocity[j], -1)).values()))
    checks.check("material-label change does not remove actual ambient stretching",
                 any(s.simplify(value) != 0 for value in f.derivative(f.mul(zeta, tag), 2).values()))
    checks.check("periodic inverse curl uses the actual base field",
                 inverse_curl(zeta) == base)

    positive_excess = f.add(zeta, {f.ZERO: s.Integer(3)})
    source = f.mul(zeta, tag)
    source_velocity = inverse_curl(source)
    reciprocal_sum = [f.add(f.mul(positive_excess, source_velocity[j]),
                            f.mul(source, base[j])).get(mean_wave, 0) for j in range(2)]
    checks.check("complete periodic inverse-curl reciprocity fixes the translation stretching row",
                 all(s.simplify(value) == 0 for value in reciprocal_sum))
    b = f.mul(carrier, base[0])
    db = gradient_green(b)
    jv = (f.scale(base[1], -1), base[0])
    potential_sum = [f.add(f.mul(positive_excess, db[j]), f.mul(jv[j], b))
                     .get(mean_wave, 0) for j in range(2)]
    checks.check("complete periodic gradient reciprocity fixes the transverse mean-current row",
                 all(s.simplify(value) == 0 for value in potential_sum)
                 and f.mul(jv[1], b).get(mean_wave, 0) != 0)

    gamma, area = s.symbols("Gamma Ac", positive=True)
    m, mz, vrchi, rchi_dot, rz = s.symbols("m Vb ZVchi Rchi_dot Rb")
    # A single component of J is represented by the commuting symbol j.
    j = s.Symbol("j")
    adot = m-s.I*k*vrchi/gamma-s.I*k*j*mz/gamma
    mdot = -s.I*k*mz/area
    # d_t integral zeta r chi = integral zeta r b+integral zeta v chi.
    corrected = adot+s.I*k*(rz+vrchi)/gamma-area*j*mdot/gamma
    checks.check("finite-k core translation keeps its actual ambient production current",
                 s.expand(corrected-m-s.I*k*rz/gamma) == 0)

    lap, adv, op_p, op_r, op_s, op_b, op_c, op_v, op_t = s.symbols("s A P R S B C V T")
    # Exact four-block determinant with a translation core L_c=0;
    # no pressure/reaction entry is prescribed from the desired output.
    matrix = s.Matrix([[lap, 0, 0, s.I*k*op_v],
                       [-op_t, lap, -s.I*k*op_b, -s.I*k*op_c],
                       [0, 0, lap+adv, -1],
                       [0, s.I*k*op_p, -k*k*op_r, lap+adv-k*k*op_s]])
    ra = 1/(lap+adv)
    wop = lap+adv-k*k*(op_r*ra+op_s)-k*k*op_p*(op_b*ra+op_c)/lap
    expected = (lap+adv)*(lap*lap*wop+k*k*op_v*op_p*op_t)
    checks.check("full material-label Schur determinant retains all finite-wave feedback signs",
                 s.factor(matrix.det()-expected) == 0)
    row_d = s.Symbol("d")
    checks.check("the translation primitive removes its apparent transport resolvent pole",
                 s.factor(row_d*adv/(lap+adv)-row_d*(1-lap/(lap+adv))) == 0)

    # Exact compactness ladder on an arbitrary fixed transport spectral
    # atom away from zero. The infinite-dimensional upgrade is the
    # finite-rank approximation proof, not a finite spectrum computation.
    eps, frequency, z = s.symbols("epsilon frequency z", positive=True)
    checks.check("each nonzero transport spectral atom has the required strong resolvent limit",
                 s.limit(eps/(eps*z+s.I*frequency), eps, 0) == 0)
    checks.check("a zero-frequency atom would change the limiting claim and is exposed",
                 s.limit(eps/(eps*z), eps, 0) == 1/z)
    cv = s.Symbol("Cv", positive=True)
    initial = s.Symbol("V0")
    mean_resolvent = lap*initial/(lap*lap+k*k*cv)
    checks.check("physical acoustic rescaling retains the mass-one common-velocity preparation",
                 s.simplify((k*mean_resolvent.subs(lap, k*z))
                            -z*initial/(z*z+cv)) == 0)
    print("Scope: exact finite-k Euler/label and current/Schur identities; no inverse-Laplace tail bound is inferred.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
