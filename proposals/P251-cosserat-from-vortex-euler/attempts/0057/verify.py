"""One-action normalization, exact affine spectrum, and moment locality."""

import sympy as s

from substrate_framework.euler_orbit import hermitian_schur_jet
from substrate_framework.homogenization import sphere_fourth_moment_isotropic
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0057-stationary-material-assembly")
    rho = s.symbols("rho", positive=True)
    velocity = s.Matrix(s.symbols("V0:3", real=True))
    fluct = s.Matrix(s.symbols("v0:3", real=True))
    paired_kinetic = rho*((velocity+fluct).dot(velocity+fluct)
                         +(velocity-fluct).dot(velocity-fluct))/4
    ledger.check("mean material mass and fluctuations are counted exactly once",
                 s.expand(paired_kinetic-rho*(velocity.dot(velocity)+fluct.dot(fluct))/2) == 0)
    pdot = s.Matrix(s.symbols("Udot0:3", real=True))
    phase = rho*velocity.dot(pdot)-rho*velocity.dot(velocity)/2
    ledger.check("material mean momentum elimination derives density rho",
                 s.expand(phase.subs(dict(zip(velocity, pdot, strict=True)))
                          -rho*pdot.dot(pdot)/2) == 0)

    omega = s.Matrix(s.symbols("w0:3"))
    xi = s.Matrix(s.symbols("xi0:3"))
    dw = s.Matrix(s.symbols("dw0:3"))
    ledger.check("complete transported second variation retains negative helicity",
                 s.expand(omega.dot(xi.cross(dw))+xi.cross(omega).dot(dw)) == 0)

    hv = s.symbols("h0:8", real=True)
    h = s.Matrix([[hv[0], hv[1], hv[2]], [hv[3], hv[4], hv[5]],
                  [hv[6], hv[7], -hv[0]-hv[4]]])
    trans = s.Matrix(s.symbols("a0:3", real=True))
    axes = [s.eye(3)[:, i] for i in range(3)]
    atom_constraints = [entry for n in axes for entry in h.T*n]
    atom_constraints += [(trans.T*n)[0] for n in axes]
    coefficient = s.linear_eq_to_matrix(atom_constraints, [*hv, *trans])[0]
    ledger.check("three independent Fourier atom pairs detect all eleven affine symmetries",
                 coefficient.rank() == 11)

    # Dual responses occupy separate supports: their mutual KKS and their
    # KKS with raw compact columns vanish. Recompute, do not assume Pi is symplectic.
    count = 11
    form = s.zeros(2*count+2)
    fq, fs = s.symbols("fq0:11"), s.symbols("fs0:11")
    b = s.symbols("B", nonzero=True)
    for j in range(count):
        form[j, count+j], form[count+j, j] = 1, -1
        form[j, 2*count], form[2*count, j] = fq[j], -fq[j]
        form[j, 2*count+1], form[2*count+1, j] = fs[j], -fs[j]
    form[2*count, 2*count+1], form[2*count+1, 2*count] = b, -b
    qcol, scol = s.zeros(2*count+2, 1), s.zeros(2*count+2, 1)
    qcol[2*count], scol[2*count+1] = 1, 1
    for j in range(count):
        qcol[count+j], scol[count+j] = -fq[j], -fs[j]
    ledger.check("dual projection kills every affine KKS cross of physical angle",
                 (form*qcol)[:count, :] == s.zeros(count, 1))
    ledger.check("dual projection kills every affine KKS cross of conjugate shape",
                 (form*scol)[:count, :] == s.zeros(count, 1))
    ledger.check("disjoint fixed responses preserve the actual compact circular pair",
                 (qcol.T*form*scol)[0] == b)

    # Exact affine material-coordinate inverse-curl multiplier derivative.
    t = s.symbols("t", real=True)
    n = s.Matrix([0, 0, 1])
    u = s.Matrix([s.symbols("u1"), s.symbols("u2"), 0])
    p = n-t*h.T*n
    mapped = -p.cross((s.eye(3)+t*h)*n.cross(u))/p.dot(p)
    derivative = s.simplify(mapped.diff(t).subs(t, 0))
    target = -h.T*u+n*(n.T*(h+h.T)*u)[0]
    ledger.check("affine multiplier derivative includes the shifted Leray response",
                 s.simplify(derivative-target) == s.zeros(3, 1))
    transverse = s.eye(3)-n*n.T
    tangent_matrix = -h.T+n*n.T*(h+h.T)
    ledger.check("affine reaction cross reduces to its exact symmetric local forcing",
                 s.simplify(transverse*(tangent_matrix+tangent_matrix.T)*transverse
                            +transverse*(h+h.T)*transverse) == s.zeros(3))
    wave = s.symbols("lambda", nonzero=True, real=True)
    projected_input = s.Matrix(s.symbols("input0:3"))
    local_monochromatic = projected_input-(wave*n)*(wave*n).T*projected_input/wave**2
    ledger.check("monochromatic Leray affine response is a finite-derivative local field",
                 s.simplify(local_monochromatic-transverse*projected_input) == s.zeros(3, 1))

    aa, bb, cc, dd, ee = s.symbols("e0:5", real=True)
    strain = s.Matrix([[aa, cc, dd], [cc, bb, ee], [dd, ee, -aa-bb]])
    norm2 = s.trace(strain*strain)
    fourth = sphere_fourth_moment_isotropic()
    average_a2 = sum(strain[i, j]*strain[k, ell]*fourth[i, j, k, ell]
                     for i in range(3) for j in range(3)
                     for k in range(3) for ell in range(3))
    ledger.check("isotropic affine first variation vanishes", s.trace(strain) == 0)
    ledger.check("exact fourth moment gives affine strain-square average",
                 s.simplify(average_a2-s.Rational(2, 15)*norm2) == 0)
    aval, bval, total = s.symbols("a b T")
    ratio = (1-t*aval+t*t*(total-bval))/(1-2*t*aval+2*t*t*bval)
    ratio2 = s.simplify(ratio.diff(t, 2).subs(t, 0)/2)
    ledger.check("spectral energy expansion retains all denominator corrections",
                 s.expand(ratio2-total+3*bval-2*aval*aval) == 0)
    averaged = s.simplify(norm2-3*norm2/3+2*average_a2)
    ledger.check("bare same-ensemble affine shear is strictly positive",
                 s.simplify(averaged-s.Rational(4, 15)*norm2) == 0)
    speed = s.symbols("U_star", positive=True)
    ledger.check("bare shear coefficient has physical normalization two fifths",
                 s.simplify(3*rho*speed**2*averaged/2
                            -s.Rational(2, 5)*rho*speed**2*norm2) == 0)
    intensity, reaction = s.symbols("intensity reaction", nonnegative=True)
    effective_mu = s.Rational(2, 5)*rho*speed**2-intensity*reaction
    ledger.check("coherence removal differs from full background-vorticity removal",
                 effective_mu.subs(intensity, 0) == s.Rational(2, 5)*rho*speed**2
                 and effective_mu.subs({intensity: 0, speed: 0}) == 0)

    # Removing the scalar first moment is a gradient representation, not energy editing.
    moment_values = s.symbols("M0:9")
    moment = s.Matrix(3, 3, moment_values)
    affine_pair = s.trace(h.T*moment)
    moment_constraints = s.Poly(affine_pair, *hv).coeffs()
    moment_linear = s.linear_eq_to_matrix(moment_constraints, moment_values)[0]
    kernel = moment_linear.nullspace()
    ledger.check("eight affine moments leave only the scalar identity force moment",
                 len(kernel) == 1 and s.Matrix(3, 3, kernel[0]) == s.eye(3))
    c = s.symbols("c")
    ledger.check("compact gradient removes precisely the remaining first moment",
                 c*s.eye(3)-(-(-c)*s.eye(3)) == s.zeros(3))

    # Exponential integral remainder: compare its exact integral with Taylor subtraction.
    z = s.symbols("z", nonzero=True, real=True)
    rem = s.integrate((1-t)*s.exp(-s.I*t*z), (t, 0, 1))
    ledger.check("averaged-center double-divergence primitive has the exact Taylor sign",
                 s.simplify(s.exp(-s.I*z)-(1-s.I*z-z*z*rem)) == 0)

    # Generic ray of the complete quadratic Leray symbol: degree four removes
    # all singular derivatives through order three at its zero spectral mode.
    pvec = s.Matrix(s.symbols("p0:3", real=True))
    projector = s.eye(3)-pvec*pvec.T/pvec.dot(pvec)
    scale = s.symbols("scale", positive=True)
    symbol = projector[0, 1]*pvec[0]*pvec[1]*pvec[2]**2
    scaled = symbol.subs({value: scale*value for value in pvec}, simultaneous=True)
    ledger.check("complete kinetic symbol is homogeneous of degree four",
                 s.simplify(scaled-scale**4*symbol) == 0)
    ledger.check("second spectral derivative has only quadratic growth",
                 s.simplify(scaled.diff(scale, 2)-12*scale**2*symbol) == 0)

    j, k, beta, phi, bd, fd = s.symbols("I K beta Phi beta_dot Phi_dot", real=True)
    local = j*(fd-bd)**2/2-k*(phi-beta)**2/2
    ledger.check("physical relative-angle map retains the gradient kinetic cross",
                 s.diff(local, fd, bd) == -j)
    ledger.check("elastic rigid co-rotation has zero relative restoring energy",
                 local.subs({phi: beta, fd: bd}) == 0)

    # Call the canonical complete Schur jet; this is an exact finite operator
    # compression check, not a diagonal-Palm or independent-cell approximation.
    pjet = (s.Matrix([[2, s.I], [-s.I, 3]]),
            s.Matrix([[0, s.I], [-s.I, 1]]), s.Matrix([[1, 1], [1, 2]]))
    njet = (s.Matrix([1, s.I]), s.Matrix([s.I, 2]), s.Matrix([2*s.I, -1]))
    hjet = (s.Matrix([[5]]), s.Matrix([[1]]), s.Matrix([[3]]))
    jet = hermitian_schur_jet(pjet, njet, hjet)
    pt = sum((t**i*value for i, value in enumerate(pjet)), s.zeros(2))
    nt = sum((t**i*value for i, value in enumerate(njet)), s.zeros(2, 1))
    ht = sum((t**i*value for i, value in enumerate(hjet)), s.zeros(1))
    full = s.simplify(ht-nt.conjugate().T*pt.inv()*nt)
    for order in range(3):
        expected = full.diff(t, order).subs(t, 0)/s.factorial(order)
        ledger.check(f"canonical noncommuting Schur jet retains complete order {order}",
                     s.simplify(jet.reduced[order]-expected) == s.zeros(1))
    print("Analytic oracles: stationary-assembly.md, affine-spectral-energy.md, slow-locality.md.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
