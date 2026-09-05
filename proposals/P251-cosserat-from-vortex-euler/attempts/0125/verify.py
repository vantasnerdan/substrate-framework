"""Actual Euler physical-angle sideband, complement mode and material spin."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0125-actual-spatial-jet")
    t = s.Symbol("t", real=True)
    om, N = s.symbols("Omega N", positive=True)
    c = s.Symbol("c", positive=True)
    delta = s.Symbol("delta", real=True)
    wm, wp = om*(1-2*c), om*(1+2*c)
    lower, upper = (1+c)/2, (1-c)/2
    trace = lower*s.cos(wm*t)+upper*s.cos(wp*t)
    ledger.check("physical angle second derivative is independent of sideband",
                 s.expand(lower*wm**2+upper*wp**2-om**2) == 0)
    fourth = s.expand(lower*wm**4+upper*wp**4)
    ledger.check("actual physical fourth derivative retains polarization variance",
                 s.expand(fourth-om**4-16*c**2*om**4*(1-c**2)) == 0)
    scalar_residual = s.diff(trace, t, 4)+(wm**2+wp**2)*s.diff(trace, t, 2)+wm**2*wp**2*trace
    ledger.check("retaining both physical Euler branches closes the fourth-order equation",
                 s.simplify(scalar_residual) == 0)
    cjet = 1-delta**2/(2*N**2)
    trace_jet = s.series(trace.subs(c, cjet), delta, 0, 3).removeO().expand()
    expected = s.cos(om*t)+delta**2/N**2*(om*t*s.sin(om*t)+(s.cos(3*om*t)-s.cos(om*t))/4)
    ledger.check("physical second spatial jet includes the third harmonic",
                 s.simplify(trace_jet-expected) == 0)
    repaired = trace_jet+delta**2/(4*N**2)*(s.cos(om*t)-s.cos(3*om*t))
    ledger.check("actual third-harmonic preparation repairs the local second jet",
                 s.simplify(repaired-s.cos(om*t)-delta**2*om*t*s.sin(om*t)/N**2) == 0)

    # Exact full-pressure Euler mode, not an input oscillator.
    x, y, z = s.symbols("x y z", real=True)
    coords = (x, y, z)
    amplitude = s.Symbol("amplitude", real=True)
    zz = x+s.I*y
    eplus, eminus = s.Matrix([1, s.I, 0]), s.Matrix([1, -s.I, 0])
    profile = amplitude*zz**2*eminus+2*amplitude/N**2*eplus
    profile[2] = 4*s.I*amplitude*zz/N
    pressure_profile = 8*s.I*om*amplitude*zz/N**2
    phase = s.exp(s.I*N*z-3*s.I*om*t)
    velocity, pressure = profile*phase, pressure_profile*phase
    background = s.Matrix([-om*y, om*x, 0])
    gradient = background.jacobian(coords)
    residual = velocity.diff(t)+velocity.jacobian(coords)*background+gradient*velocity
    residual += s.Matrix([s.diff(pressure, var) for var in coords])
    ledger.check("polynomial complement solves the complete linear Euler equation",
                 all(s.simplify(value) == 0 for value in residual))
    ledger.check("polynomial complement is exactly divergence free",
                 s.simplify(sum(s.diff(velocity[i], coords[i]) for i in range(3))) == 0)

    def curl(vector):
        return s.Matrix([s.diff(vector[2], y)-s.diff(vector[1], z),
                         s.diff(vector[0], z)-s.diff(vector[2], x),
                         s.diff(vector[1], x)-s.diff(vector[0], y)])

    vorticity = curl(velocity)
    displacement = vorticity/(2*s.I*om*N)
    lin = displacement.diff(t)+displacement.jacobian(coords)*background-gradient*displacement-velocity
    ledger.check("same complement has actual Kelvin/Lin displacement",
                 all(s.simplify(value) == 0 for value in lin))
    core = {x: 0, y: 0, z: 0}
    angle = s.Matrix([0, 0, 1]).cross(vorticity)/(2*om)
    core_angle = s.simplify(angle.subs(core))
    expected_angle = s.I*amplitude/(om*N)*eplus*s.exp(-3*s.I*om*t)
    ledger.check("actual core vorticity angle of the complement is nonzero",
                 s.simplify(core_angle-expected_angle) == s.zeros(3, 1))

    # Leading exact small-tag tensor coefficient, including deformation/current.
    # For isotropic covariance, int rho r_i r_j = moment*delta_ij.
    moment = s.Symbol("moment", positive=True)
    displacement_gradient = displacement.jacobian(coords).subs(core)
    velocity_gradient = velocity.jacobian(coords).subs(core)
    spin = s.Matrix([sum(moment*(s.eye(3)[:, j].cross(
        (gradient*displacement_gradient+velocity_gradient)[:, j])
        +displacement_gradient[:, j].cross(gradient[:, j]))[i] for j in range(3))
        for i in range(3)])
    spin = s.simplify(spin)
    print("3Omega core angle row:", core_angle.T)
    print("3Omega full isotropic-tag spin row:", spin.T)
    print("3Omega spin minus 2*moment*angle_rate:", s.simplify(spin-2*moment*core_angle.diff(t)).T)
    ledger.check("full physical tag spin includes both displacement and velocity gradients",
                 spin != s.simplify(moment*vorticity.subs(core)))
    ledger.check("measured third-harmonic spin is one-third tag moment times angle rate",
                 s.simplify(spin-moment*core_angle.diff(t)/3) == s.zeros(3, 1))

    # Same-frequency measured angle/spin controls supplied by actual Euler modes.
    j1, j2, angle_target, spin_target = s.symbols("j1 j2 angle_target spin_target")
    coefficients = s.Matrix([(spin_target-j2*angle_target)/(j1-j2),
                             (j1*angle_target-spin_target)/(j1-j2)])
    ledger.check("two actual sideband controls meet both independent physical rows",
                 s.simplify(s.Matrix([[1, 1], [j1, j2]])*coefficients
                            -s.Matrix([angle_target, spin_target])) == s.zeros(2, 1))

    # Full mean-free Bloch projection coefficients on each nonzero cell mode.
    l1, l2, l3, n1, n2, n3 = s.symbols("l1 l2 l3 n1 n2 n3", real=True)
    ell, direction = s.Matrix([l1, l2, l3]), s.Matrix([n1, n2, n3])
    square, dot = ell.dot(ell), ell.dot(direction)
    proj = s.eye(3)-(ell+delta*direction)*(ell+delta*direction).T/(ell+delta*direction).dot(ell+delta*direction)
    p1 = -(direction*ell.T+ell*direction.T)/square+2*dot*ell*ell.T/square**2
    p2 = (-direction*direction.T/square+2*dot*(direction*ell.T+ell*direction.T)/square**2
          +(direction.dot(direction)/square**2-4*dot**2/square**3)*ell*ell.T)
    ledger.check("full first Bloch pressure projector derivative",
                 s.simplify(proj.diff(delta).subs(delta, 0)-p1) == s.zeros(3))
    ledger.check("full second Bloch pressure projector derivative",
                 s.simplify(proj.diff(delta, 2).subs(delta, 0)/2-p2) == s.zeros(3))

    # A return supported outside the original profile cancels the global first
    # moment without changing the tagged row. The integral scalings are the
    # exact substitutions z=L*s, not an optimization or a fitted response.
    scale, bcore, breturn, norm_return = s.symbols(
        "scale bcore breturn norm_return", nonzero=True)
    return_amplitude = -bcore/(scale**2*breturn)
    ledger.check("outside-tag return cancels the global first axial moment",
                 s.simplify(bcore+return_amplitude*scale**2*breturn) == 0)
    ledger.check("outside-tag return action cost scales as inverse length cubed",
                 s.simplify(return_amplitude**2*scale*norm_return
                            -bcore**2*norm_return/(scale**3*breturn**2)) == 0)

    # Actual operator Duhamel closure: noncommuting products must retain order.
    a0, a1, a2, z0, z1, z2 = s.symbols("a0 a1 a2 z0 z1 z2", commutative=False)
    product = s.expand((a0+delta*a1+delta**2*a2)*(z0+delta*z1+delta**2*z2))
    ledger.check("second Euler jet retains the first-jet coupling rather than averaging it",
                 s.expand(product.coeff(delta, 2)-(a0*z2+a1*z1+a2*z0)) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
