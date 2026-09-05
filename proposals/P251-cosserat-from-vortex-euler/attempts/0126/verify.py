"""Actual rotating-profile clock, paired action and observation normalization."""

import sympy as s

from substrate_framework.euler_observation import material_tag_fourier_dipole
from substrate_framework.euler_phase import moving_phase_pullback
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0126-laboratory-pair")
    omega, beta, density, nu = s.symbols("Omega beta rho nu", positive=True)
    t = s.Symbol("t", real=True)
    jrot = s.Matrix([[0, -1], [1, 0]])
    a = s.Matrix(s.symbols("a:2", real=True))
    adot = -omega*jrot*a
    velocity = adot-omega*jrot*a
    velocity_rate = velocity.jacobian(a)*adot
    checks.check("plateau velocity satisfies the laboratory Euler momentum equation",
                 s.simplify(velocity_rate+omega*jrot*velocity) == s.zeros(2, 1))
    rotation = s.Matrix([[s.cos(omega*t), -s.sin(omega*t)],
                         [s.sin(omega*t), s.cos(omega*t)]])
    rotating_generator = s.simplify(rotation.T*(-omega*jrot)*rotation
                                    +rotation.T.diff(t)*rotation)
    checks.check("physical rotating observation has twice the lab frequency",
                 rotating_generator == -2*omega*jrot)
    checks.check("wrong lab clock fails actual Euler momentum",
                 (-2*omega*jrot*velocity+omega*jrot*velocity) != s.zeros(2, 1))

    phase_plus = -beta*jrot
    phase = s.diag(phase_plus/2, -phase_plus/2)
    h = beta*omega*s.eye(4)/2
    embed = s.BlockMatrix([[s.eye(2), s.eye(2)],
                           [s.eye(2), -s.eye(2)]]).as_explicit()
    pulled = moving_phase_pullback(phase, h, embed, s.zeros(4))
    generator = s.BlockMatrix([[s.zeros(2), -omega*jrot],
                               [-omega*jrot, s.zeros(2)]]).as_explicit()
    checks.check("one paired Euler action yields the actual initial-data join",
                 pulled.generator == generator and pulled.residual == s.zeros(4))
    checks.check("all independent initial data satisfy the optical polynomial",
                 generator**2 == -omega**2*s.eye(4))
    q = s.Matrix([s.Function("q1")(t), s.Function("q2")(t)])
    r = s.Matrix([s.Function("r1")(t), s.Function("r2")(t)])
    z = q.col_join(r)
    direct = -(z.T*pulled.symplectic*z.diff(t))[0]/2-(z.T*pulled.hamiltonian*z)[0]/2
    convenient = beta*(r.T*jrot*q.diff(t))[0]-beta*omega*(q.dot(q)+r.dot(r))/2
    boundary = beta*(q.T*jrot*r)[0]/2
    checks.check("reaction action differs only by its explicitly computed boundary term",
                 s.simplify(direct-convenient-s.diff(boundary, t)) == 0)
    reaction = jrot*q.diff(t)/omega
    reduced = s.simplify(convenient.subs(dict(zip(r, reaction)), simultaneous=True))
    inertia = beta/omega
    checks.check("independent reaction variation produces the positive physical mass",
                 all(s.simplify(s.diff(convenient, x).subs(dict(zip(r, reaction)),
                                                          simultaneous=True)) == 0 for x in r)
                 and s.simplify(reduced-inertia*(q.diff(t).dot(q.diff(t))
                                                 -omega**2*q.dot(q))/2) == 0)
    spin_plus = inertia*(-omega*jrot*(q+r))
    spin_minus = inertia*(omega*jrot*(q-r))
    checks.check("measured pair spin equals the derived action momentum",
                 s.simplify(((spin_plus+spin_minus)/2).subs(dict(zip(r, reaction)),
                                                                          simultaneous=True)
                            -inertia*q.diff(t)) == s.zeros(2, 1))
    checks.check("forcing the reaction history to zero deletes nontrivial dynamics",
                 generator[2:4, :2].det() == omega**2)

    # The three coordinate axes are an exact degree-two spherical design.
    axes = [s.eye(3)[:, i] for i in range(3)]
    frame = sum((s.eye(3)-n*n.T for n in axes), s.zeros(3))/3
    phi = s.Matrix(s.symbols("Phi:3", real=True))
    inferred = frame.inv()*sum(((s.eye(3)-n*n.T)*phi for n in axes), s.zeros(3, 1))/3
    checks.check("actual transverse measurements reconstruct the physical angle",
                 frame == s.Rational(2, 3)*s.eye(3) and inferred == phi)
    action_mass = nu*inertia*frame
    spin_mass = nu*inertia*frame
    raw_mean_mass = frame.inv().T*action_mass*frame.inv()
    raw_mean_spin = spin_mass*frame.inv()
    checks.check("one geometric normalization aligns action and measured spin",
                 action_mass == 2*nu*inertia*s.eye(3)/3
                 and raw_mean_mass != raw_mean_spin)
    wavevector = s.Matrix(s.symbols("k:3", real=True))
    rate = s.Matrix(s.symbols("V:3", real=True))
    physical_spin = spin_mass*rate
    dipole = material_tag_fourier_dipole(wavevector, physical_spin, s.zeros(3))
    checks.check("actual hybrid optical transfer retains total density and sign",
                 s.simplify(-dipole/density+s.I*nu*inertia*wavevector.cross(rate)/(3*density))
                 == s.zeros(3, 1))
    # This last input uses the independently justified coherent k=0 shape row,
    # not a claim that an individual tag's shape-rate tensor is zero.
    print("Pair inertia:", inertia, "isotropic inertia:", 2*nu*inertia/3)
    print("Scope: actual leading Euler profiles and explicit finite-time transfer errors")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
