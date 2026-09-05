"""Same-cell body metric, first shape connection and covariant field map."""

import sys
from pathlib import Path

import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "0040"))
from fourier_orbit import (  # noqa: E402
    actual_tube, core_cage_generators, cross, derivative, mul,
    orbit_matrices, scale, trig,
)
from substrate_framework.verification import CheckLedger  # noqa: E402


def main():
    ledger = CheckLedger("P251-0042-body-map")
    rho, ell, b = s.symbols("rho ell b", positive=True)
    x, y, z = s.symbols("x y z", real=True)
    body_position = ell*s.Matrix([x, y, z])
    fluid_velocity = s.Matrix([-b*s.sin(y), 2*b*s.sin(x), 2*b*s.cos(x)+b*s.cos(y)])
    spin_generator = s.Matrix([-body_position[1], body_position[0], 0])
    omega = s.symbols("Omega", real=True)
    kinetic_increment = rho*((fluid_velocity+omega*spin_generator).dot(
        fluid_velocity+omega*spin_generator)-fluid_velocity.dot(fluid_velocity))/2
    angular_density = rho*body_position.cross(fluid_velocity)[2]
    inertia_density = rho*spin_generator.dot(spin_generator)
    ledger.check("body kinetic decomposition retains the full angular-momentum connection",
                 s.expand(kinetic_increment-omega*angular_density-omega**2*inertia_density/2) == 0)

    def cell_mean(expr):
        return s.simplify(s.integrate(s.integrate(s.integrate(expr,
                             (z, -s.pi, s.pi)), (y, -s.pi, s.pi)),
                             (x, -s.pi, s.pi))/(2*s.pi)**3)

    body_inertia = cell_mean(inertia_density)
    background_momentum = cell_mean(angular_density)
    ledger.check("same cube has computed locked inertia 2*rho*pi^2*ell^2/3",
                 body_inertia == 2*rho*s.pi**2*ell**2/3)
    ledger.check("background angular momentum is retained before removing its total derivative",
                 background_momentum == 3*rho*b*ell)
    ledger.check("stationary field has zero mean linear velocity",
                 fluid_velocity.applyfunc(cell_mean) == s.zeros(3, 1))

    angle, _ = core_cage_generators(1)
    f = mul(trig(0, kind="sin"), trig(1, kind="sin"))
    partner = (mul(derivative(f, 0), trig(2, kind="sin")),
               mul(derivative(f, 1), trig(2, kind="sin")),
               scale(mul(f, trig(2)), -2))
    generators = (angle, partner)
    tangents, _, _ = orbit_matrices(actual_tube(2*b, b), generators)
    for index in range(2):
        ledger.check(f"shape {index} has zero axial average at every cross-section point",
                     all(all(wave[2] != 0 for wave in component)
                         for component in generators[index]))
        ledger.check(f"physical velocity tangent {index} has zero axial average",
                     all(all(wave[2] != 0 for wave in component)
                         for component in tangents[index]))
    # Moment weights x,y and the base field are independent of z. The preceding
    # exact zero axial means therefore also control moving-parcel first variations.
    rotation_shape = cross(generators[0], generators[1])[2]
    ledger.check("axial shape-connection curvature has zero mean",
                 rotation_shape.get((0, 0, 0), 0) == 0)
    qshape = ell*s.Matrix([-s.sin(y)*s.cos(z), s.sin(x)*s.cos(z), 0])
    sshape = ell*s.Matrix([s.cos(x)*s.sin(y)*s.sin(z),
                           s.sin(x)*s.cos(y)*s.sin(z),
                           -2*s.sin(x)*s.sin(y)*s.cos(z)])
    for index, generator in enumerate((qshape, sshape)):
        ledger.check(f"body-shape kinetic cross metric {index} vanishes",
                     cell_mean(rho*spin_generator.dot(generator)) == 0)
        ledger.check(f"first material axial-inertia variation {index} vanishes",
                     cell_mean(2*rho*(body_position[0]*generator[0]
                                     +body_position[1]*generator[1])) == 0)

    # Routh/body extension: exact algebra once the transported body coordinate is licensed.
    iq = 6*rho*ell**2/25
    kq = 7*rho*b**2/48
    beta_dot, psi_dot, psi, beta = s.symbols("beta_dot Psi_dot Psi beta", real=True)
    qdot, q = s.symbols("q_dot q", real=True)
    coefficient = s.factor(1+iq/body_inertia)
    kinetic = body_inertia*(beta_dot+qdot)**2/2+iq*qdot**2/2
    mapped = s.expand(kinetic.subs(qdot, (psi_dot-beta_dot)/coefficient))
    jpsi = s.factor(body_inertia**2/(body_inertia+iq))
    jbeta = s.factor(body_inertia*iq/(body_inertia+iq))
    kpsi = s.factor(kq/coefficient**2)
    ledger.check("computed coefficient makes the complete kinetic map diagonal",
                 s.simplify(mapped-(jpsi*psi_dot**2+jbeta*beta_dot**2)/2) == 0)
    ledger.check("same coordinate map transports the actual restoring potential",
                 s.simplify((kq*q*q/2).subs(q, (psi-beta)/coefficient)
                            -kpsi*(psi-beta)**2/2) == 0)
    common = s.symbols("common_rotation", real=True)
    fieldmap = beta+coefficient*q
    ledger.check("mapped microrotation has unit common-rotation weight",
                 s.simplify(fieldmap.subs(beta, beta+common)-fieldmap-common) == 0)
    ledger.check("fixed-cage frequency includes the body's computed inertia",
                 s.simplify(kpsi/jpsi-kq/(body_inertia+iq)) == 0)
    ledger.check("fixed-cage and fixed-body frequencies are physically distinct",
                 s.simplify(kpsi/jpsi-kq/iq) != 0)

    # Expose the RVE/material-parcel difference using its actual Euler moment flux.
    qxy_dot = cell_mean(rho*(body_position[0]*fluid_velocity[1]
                             +body_position[1]*fluid_velocity[0]))
    ledger.check("an advected cube immediately changes its off-diagonal moment",
                 qxy_dot == rho*ell*b)
    ledger.check("the fixed periodic face is not a material surface",
                 fluid_velocity[0].subs(x, s.pi) == -b*s.sin(y))
    # Macro expansion about cell centers counts its affine rotation only once.
    velocity = s.Matrix(s.symbols("V0:3"))
    affine = s.Matrix(3, 3, lambda i, j: s.Symbol(f"A{i}{j}", real=True))
    affine_velocity = velocity+affine*body_position
    affine_energy = cell_mean(rho*affine_velocity.dot(affine_velocity)/2)
    expected = rho*velocity.dot(velocity)/2+rho*s.pi**2*ell**2*sum(
        value**2 for value in affine)/6
    ledger.check("cell-center affine kinetic expansion computes the finite-cell correction",
                 s.simplify(affine_energy-expected) == 0)
    ledger.check("its antisymmetric part is exactly the body rotation energy",
                 s.simplify((rho*s.pi**2*ell**2*sum(value**2 for value in affine)/6)
                 .subs(dict(zip(affine, s.Matrix([[0, -omega, 0], [omega, 0, 0], [0, 0, 0]]))))
                            -body_inertia*omega**2/2) == 0)
    print("I_body per volume:", body_inertia)
    print("L0 per volume:", background_momentum)
    print("covariant coefficient A:", coefficient)
    print("J_Psi per volume:", jpsi)
    print("J_beta per volume:", jbeta)
    print("K_Psi per volume:", kpsi)
    print("fixed-cage omega_squared:", s.factor(kpsi/jpsi))
    print("Scope: exact material metric and conditional body/internal Routh field map.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
