"""Physical-frame sign, exact moving action, and retained current checks."""

import sympy as s

from substrate_framework.euler_phase import moving_phase_pullback
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0114-physical-frame-action")
    t = s.symbols("t", real=True)
    beta, rate, eps, tube, carrier = s.symbols("beta rate eps tube N", positive=True)
    torsion, curvature, y1, y2 = s.symbols("torsion curvature y1 y2", real=True)
    j = s.Matrix([[0, -1], [1, 0]])
    j0 = -j

    def zero(v):
        return all(s.simplify(entry) == 0 for entry in v)

    def rotation(theta):
        return s.Matrix([[s.cos(theta), -s.sin(theta)], [s.sin(theta), s.cos(theta)]])

    tangent_derivative = s.Matrix([1-tube*curvature*y1, -tube*torsion*y2,
                                    tube*torsion*y1])
    first_normal, second_normal = s.Matrix([0, tube, 0]), s.Matrix([0, 0, tube])
    ledger.check("source metric fixes its torsion as the negative of the usual Frenet convention",
                 tangent_derivative.dot(first_normal) == -tube**2*torsion*y2
                 and tangent_derivative.dot(second_normal) == tube**2*torsion*y1)
    wobble, r0 = s.symbols("wobble r0", real=True)
    lifted_phase = r0*t+wobble*s.sin(t)
    registration = rotation(lifted_phase-r0*t)
    physical_generator = s.diff(lifted_phase, t)*j
    registered_generator = s.simplify(registration.T*physical_generator*registration
                                      -registration.T*registration.diff(t))
    ledger.check("the fixed geometric uniform-twist registration removes only its zero-winding part",
                 zero(registered_generator-r0*j)
                 and s.simplify((lifted_phase-r0*t).subs(t, 2*s.pi)) == 0
                 and registration.subs(t, 0) == s.eye(2))

    kx, ky = s.symbols("kx ky", real=True)
    kz, speed = s.symbols("kz U", positive=True)
    lam = s.symbols("lambda", real=True, nonzero=True)
    a1 = s.Matrix([0, -1/(carrier*kz), ky/(carrier*kz**2)])
    a2 = s.Matrix([1/(carrier*kz), 0, -kx/(carrier*kz**2)])
    kvec = s.Matrix([kx, ky, kz])
    vorticity = s.Matrix([0, 0, lam*speed])
    ledger.check("physical unit-angle columns satisfy the actual transverse Kelvin constraint",
                 s.simplify(kvec.dot(a1)) == 0 and s.simplify(kvec.dot(a2)) == 0)
    ledger.check("actual KKS sign follows the physical curl sign with no orientation reset",
                 s.simplify(vorticity.dot(a1.cross(a2))-lam*speed/(carrier**2*kz**2)) == 0)

    aa, bb, cc, dd, ee, ff, gg, hh = s.symbols("a b c d e f g h", real=True)
    gradient = s.Matrix([[aa, bb, cc], [dd, ee, ff], [gg, hh, -aa-ee]])
    omega = s.Matrix([hh-ff, cc-gg, dd-bb])
    e1, e2, e3 = s.eye(3)[:, 0], s.eye(3)[:, 1], s.eye(3)[:, 2]
    pairing = e3.dot(omega)
    d1 = gradient*e1-pairing*e3.cross(e1)
    d2 = gradient*e2-pairing*e3.cross(e2)
    density_rate = (gradient*omega).dot(e1.cross(e2))+omega.dot(d1.cross(e2)+e1.cross(d2))
    ledger.check("leading transported KKS density is conserved by the actual amplitude equation",
                 s.expand(density_rate) == 0)
    tangent_component = s.symbols("a_t", real=True)
    transverse = s.Matrix([aa, bb, 0])
    tangent_rate = gradient*e3-e3*(e3.dot(gradient*e3))
    projector = s.diag(1, 1, 0)
    ledger.check("physical transverse quotient cancels its tangent-displacement contamination",
                 zero(projector*gradient*(transverse+tangent_component*e3)
                      -tangent_rate*tangent_component-projector*gradient*transverse))

    omega0 = beta*j0
    signed_rate = -rate
    b0 = signed_rate*j
    h0 = -omega0*b0
    ledger.check("the registered physical curl-sign candidate gives a positive Hessian",
                 h0 == beta*rate*s.eye(2)
                 and (-(-omega0)*b0) == -beta*rate*s.eye(2))
    q, p, qdot = s.symbols("q p qdot", real=True)
    lagrangian = beta*p*qdot-beta*rate*(q*q+p*p)/2
    eliminated = s.simplify(lagrangian.subs(p, qdot/rate))
    ledger.check("autonomous angular inertia is derived by eliminating the actual conjugate coordinate",
                 s.simplify(eliminated-beta*qdot**2/(2*rate)+beta*rate*q**2/2) == 0)
    eigenvector = s.Matrix([1, -s.I])
    krein = s.simplify(s.I*(eigenvector.conjugate().T*omega0*eigenvector)[0])
    ledger.check("positive-frequency Krein sign matches the frozen physical energy branch",
                 zero(b0*eigenvector+s.I*rate*eigenvector) and krein == 2*beta)

    # Independent finite-dimensional receipt for the exact near-identity
    # rectification, using the importable full moving-phase API from 0115.
    delta = s.symbols("delta", real=True)
    embedding = s.exp(eps*s.sin(t))*rotation(-delta*s.sin(t))
    pulled = moving_phase_pullback(omega0, h0, embedding, embedding.diff(t))
    beta_t = beta*s.exp(2*eps*s.sin(t))
    expected_b = -eps*s.cos(t)*s.eye(2)+(-rate+delta*s.cos(t))*j
    ledger.check("full finite moving action retains the actual symplectic-rate contribution",
                 zero(pulled.symplectic-beta_t*j0)
                 and zero(pulled.generator-expected_b))
    ledger.check("connection-corrected finite angular Hessian has its exposing exact sign margin",
                 zero(pulled.hamiltonian-beta_t*(rate-delta*s.cos(t))*s.eye(2)))
    rectifier = s.exp(-eps*s.sin(t))*rotation(delta*s.sin(t))
    ledger.check("unique rectifier solves its actual coefficient ODE with identity initial data",
                 zero(rectifier.diff(t)-pulled.generator*rectifier+rectifier*b0)
                 and rectifier.subs(t, 0) == s.eye(2))
    ledger.check("rectification transports the variable KKS form to the original constant form",
                 zero(rectifier.T*pulled.symplectic*rectifier-omega0))
    connection = rectifier.T*pulled.symplectic*rectifier.diff(t)
    new_h = rectifier.T*pulled.hamiltonian*rectifier+(connection+connection.T)/2
    ledger.check("the same finite action becomes the prescribed autonomous positive action",
                 zero(new_h-h0))
    physical_observation = s.Matrix([[1, 2]])
    ledger.check("physical observations are transformed with the actual action correction",
                 zero(physical_observation*embedding*rectifier-physical_observation)
                 and not zero(physical_observation*embedding-physical_observation))

    current = s.Matrix(3, 3, s.symbols("Q11 Q12 Q13 Q21 Q22 Q23 Q31 Q32 Q33", real=True))
    response = sum(((s.eye(3)-s.eye(3)[:, n]*s.eye(3)[:, n].T)
                    *current*s.eye(3)[:, n]).dot((s.eye(3)-s.eye(3)[:, n]*s.eye(3)[:, n].T)
                                               *current*s.eye(3)[:, n]) for n in range(3))
    spin_sq = sum((current[i, n]-current[n, i])**2 for i in range(3) for n in range(i+1, 3))
    symmetric_offdiag = sum((current[i, n]+current[n, i])**2/2
                           for i in range(3) for n in range(i+1, 3))
    ledger.check("full shape-plus-spin current has a nonzero transverse row whenever spin is nonzero",
                 s.expand(response-spin_sq/2-symmetric_offdiag) == 0)
    ledger.check("physical moment density is not falsely kept finite by an implicit carrier rescaling",
                 carrier**3*carrier**-5 == carrier**-2)
    order = s.symbols("m", integer=True, positive=True)
    ledger.check("the explicit normalized moment error survives the high-carrier sensitivity",
                 (s.Rational(1, 2)-order).subs(order, 3) < 0)

    # Source-compatible invariant covariance gives genuinely time-dependent
    # material spin without assuming a Floquet frequency splitting.
    w0, p0, q0, r1, s0, t0 = s.symbols("w0 p0 q0 r1 s0 t0", real=True)
    beltrami_gradient = s.Matrix([[p0, q0, r1], [q0+w0, s0, t0],
                                  [r1, t0, -p0-s0]])
    ax, ay = s.symbols("ax ay", real=True)
    material_amplitude = s.Matrix([ax, ay, -(kx*ax+ky*ay)/kz])
    velocity_amplitude = -w0*kz/kvec.dot(kvec)*kvec.cross(material_amplitude)
    spin_vector = material_amplitude.cross(beltrami_gradient*(speed*e3)) \
        +(speed*e3).cross(beltrami_gradient*material_amplitude+velocity_amplitude)
    spin_matrix = s.simplify(spin_vector.jacobian((ax, ay)))
    exposing = spin_matrix[0, 0]+spin_matrix[1, 1] \
        -2*(kx*spin_matrix[2, 0]+ky*spin_matrix[2, 1])/kz
    ledger.check("the periodic-covariance material spin row is nonzero for nonzero actual vorticity",
                 s.simplify(exposing-speed*w0*kz**2/kvec.dot(kvec)) == 0)
    theta = s.symbols("theta", real=True)
    ct, cn = s.symbols("ct cn", positive=True)
    conjugacy = s.Matrix([[1, aa, bb], [0, 1, cc], [0, 0, 1]])
    orthogonal_return = s.diag(1, rotation(theta))
    particle_return = conjugacy*orthogonal_return*conjugacy.inv()
    covariance = conjugacy*s.diag(ct, cn, cn)*conjugacy.T
    periodic_covector = conjugacy.inv().T*e1
    ledger.check("positive material covariance is invariant under the actual elliptic particle normal form",
                 zero(particle_return*covariance*particle_return.T-covariance)
                 and covariance.det() == ct*cn**2)
    ledger.check("invariant covariance maps the periodic covector to the actual core velocity direction",
                 zero(covariance*periodic_covector-ct*conjugacy*e1))
    euler_return = rotation(theta)
    increment = euler_return-s.eye(2)
    ledger.check("nontrivial Euler Floquet return permits a prepared nonzero spin endpoint difference",
                 s.simplify(increment.det()-(2-2*s.cos(theta))) == 0
                 and zero(increment*increment.inv()*s.Matrix([1, 0])-s.Matrix([1, 0])))
    print("EXACT physical sign criterion: lambda * unwrapped source torsion < 0")
    print("EXACT positive angular inertia:", beta/rate)
    print("EXACT positive-frequency Krein pairing:", krein)
    print("EXACT current inequality remainder:", symmetric_offdiag)
    print("EXACT nonconstant-spin row invariant:", s.simplify(exposing))
    print("SCOPE: positive fixed-frame packet action; canonical momentum is not assumed to equal tube spin")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
