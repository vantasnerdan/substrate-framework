"""Independent material-flow check of the exact spin-correlation trace.

An explicit integrable radial shear flow checks the kinematic identity
independently of the Euler-column proof. It is NOT a replacement Beltrami
supplier. Gaussian moments here are exact probability integrals, not a
sampled tag or an approximation to the fixed C016 tag.
"""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0235-material-flow-trace")
    x, y, z, t, tau = s.symbols("x y z t tau", real=True)
    frequency, width = s.symbols("frequency width", positive=True)
    eye = s.eye(3)
    # Normal reference density exp(-y²-z²)/pi, so E[r²]=1.
    psi = -frequency*(y*y+z*z)/2
    mean = -frequency/2

    def normal_flow(time):
        return s.Matrix([y*s.cos(frequency*time)+z*s.sin(frequency*time),
                         -y*s.sin(frequency*time)+z*s.cos(frequency*time)])

    normal = normal_flow(t)
    flow = s.Matrix([x+t*psi, normal[0], normal[1]])
    center = s.Matrix([t*mean, 0, 0])
    radius = flow-center
    velocity = flow.diff(t)
    relative = velocity-center.diff(t)
    earlier = normal_flow(tau)
    lag = frequency*(t-tau)
    tangent = s.Matrix([[1, -frequency*(t-tau)*earlier[0],
                        -frequency*(t-tau)*earlier[1]],
                       [0, s.cos(lag), s.sin(lag)],
                       [0, -s.sin(lag), s.cos(lag)]])
    memory = 2*t*eye-2*tangent.applyfunc(lambda value: s.integrate(value, (tau, 0, t)))

    def expectation(value):
        result = 0
        for powers, coefficient in s.Poly(s.expand(value), x, y, z).terms():
            if any(order % 2 for order in powers):
                continue
            ix, iy, iz = powers
            result += (coefficient*s.factorial2(ix-1)*width**(ix//2)
                       *s.factorial2(iy-1)/2**(iy//2)
                       *s.factorial2(iz-1)/2**(iz//2))
        return s.simplify(s.trigsimp(result))

    # All mean/centroid correction terms vanish in the integral by centering.
    trace = 0
    for axis in range(3):
        displacement = t*eye[:, axis].cross(flow)+memory[:, axis]
        material_spin = displacement.cross(relative)+radius.cross(displacement.diff(t))
        trace += expectation(material_spin[axis])
    trace = s.simplify(s.trigsimp(trace))
    cchi = s.cos(frequency*t)
    cd = -cchi  # d=(psi-mean)chi_psi=(1-r²)chi; E[(1-r²)r²]=-1.
    inertia = width+1+frequency**2*t**2/4
    predicted = (2*inertia+2*(cchi-1)
                 +2*(t*s.diff(cd, t)-cd-1))
    print("Direct normalized material spin trace:", trace, flush=True)
    print("Defining correlation-trace defect:", s.simplify(s.trigsimp(trace-predicted)),
          flush=True)
    checks.check("direct full material-flow spin equals the complete correlation identity",
                 s.simplify(s.trigsimp(trace-predicted)) == 0)
    checks.check("the transported axial cut contributes its actual central shear variance",
                 s.simplify(expectation(radius.dot(radius))-inertia) == 0)
    checks.check("dropping either the material memory or axial cut changes the trace",
                 s.simplify(trace-2*inertia) != 0
                 and s.simplify(trace-(predicted-frequency**2*t**2/2)) != 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
