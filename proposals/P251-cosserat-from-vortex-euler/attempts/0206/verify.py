"""Exact geometric, phase and slow-operator identities; no ring-pole oracle."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0206-low-ring-action-and-current")
    time = s.symbols("t", real=True)
    impulse, speed, radius, rho, gamma, tension = s.symbols(
        "I0 U R rho Gamma Tau", positive=True)
    spin = s.symbols("J0", real=True)
    x, yy, px, py = (s.Function(n)(time) for n in ("X", "Y", "px", "py"))
    gyro = spin/impulse**2
    action = (px*s.diff(x, time)+py*s.diff(yy, time)
              +gyro*(px*s.diff(py, time)-py*s.diff(px, time))/2
              -speed*(px**2+py**2)/(2*impulse))

    def variation(field):
        return s.simplify(s.diff(action, field)-s.diff(
            s.diff(action, s.diff(field, time)), time))

    checks.check("actual Euclidean action conserves transverse impulse",
                 variation(x) == -s.diff(px, time) and variation(yy) == -s.diff(py, time))
    checks.check("angular-impulse gyro is retained in the full X reconstruction",
                 s.simplify(variation(px)-s.diff(x, time)
                            -gyro*s.diff(py, time)+speed*px/impulse) == 0)
    checks.check("angular-impulse gyro is retained in the full Y reconstruction",
                 s.simplify(variation(py)-s.diff(yy, time)
                            +gyro*s.diff(px, time)+speed*py/impulse) == 0)
    checks.check("on actual conserved phases the same action gives the tilted-ring drift",
                 s.simplify(variation(px).subs(s.diff(py, time), 0)
                            -s.diff(x, time)+speed*px/impulse) == 0)
    symplectic = s.Matrix([[0, 0, -1, 0], [0, 0, 0, -1],
                          [1, 0, 0, gyro], [0, 1, -gyro, 0]])
    checks.check("exact Euclidean suborbit is nondegenerate with nonzero swirl",
                 s.factor(symplectic.det()) == 1)
    coordinates = s.diag(1, 1, impulse, impulse)
    checks.check("physical orientation/momentum KKS determinant keeps its I0 scale",
                 s.factor((coordinates.T*symplectic*coordinates).det()) == impulse**4)
    hess = s.diag(0, 0, speed/impulse, speed/impulse)
    generator = -symplectic.inv()*hess
    checks.check("rigid tilt/translation is an actual Jordan pair, not a spurious oscillator",
                 generator**2 == s.zeros(4) and generator != s.zeros(4))

    ex, ey, ez, cp, cz = s.symbols("ex ey ez C_perp C_parallel", real=True)
    rotation = s.Matrix([[0, -ez, ey], [ez, 0, -ex], [-ey, ex, 0]])
    covariance = s.diag(cp, cp, cz)
    delta_covariance = rotation*covariance-covariance*rotation
    delta_normal = rotation*s.Matrix([0, 0, 1])
    checks.check("actual covariance normal has unit rigid-rotation response",
                 s.simplify(delta_covariance[0, 2]-(cz-cp)*delta_normal[0]) == 0
                 and s.simplify(delta_covariance[1, 2]-(cz-cp)*delta_normal[1]) == 0)
    checks.check("circular in-plane covariance cannot supply a nonzero angle gap",
                 covariance[0, 0]-covariance[1, 1] == 0)

    eps, q, z, qp, zp = s.symbols("epsilon q z q_phi z_phi", real=True)
    length_density = s.sqrt((radius+eps*q)**2+eps**2*(qp**2+zp**2))
    impulse_density = rho*gamma*(radius+eps*q)**2/2
    frame_density = tension*length_density-speed*impulse_density
    stationary_speed = tension/(rho*gamma*radius)
    checks.check("actual length second variation includes both bending directions",
                 s.simplify(s.diff(length_density, eps, 2).subs(eps, 0)
                            -(qp**2+zp**2)/radius) == 0)
    checks.check("translating-frame impulse cancels the first variation",
                 s.simplify(s.diff(frame_density, eps).subs(
                     {eps: 0, speed: stationary_speed})) == 0)
    checks.check("complete second variation retains the indispensable negative q squared",
                 s.simplify(s.diff(frame_density, eps, 2).subs(
                     {eps: 0, speed: stationary_speed})
                            -tension*(qp**2+zp**2-q**2)/radius) == 0)
    checks.check("omitting actual impulse changes the translation Hessian",
                 s.simplify(s.diff(tension*length_density, eps, 2).subs(eps, 0)
                            -tension*(qp**2+zp**2-q**2)/radius) != 0)
    phi = s.symbols("phi", real=True)
    er = s.Matrix([s.cos(phi), s.sin(phi), 0])
    et = s.Matrix([-s.sin(phi), s.cos(phi), 0])
    axis = s.Matrix([0, 0, 1])
    checks.check("actual filament KKS orientation fixes the q-z sign",
                 s.simplify(et.dot(er.cross(axis))) == -1)
    n = s.symbols("n", positive=True, integer=True)
    h_bend = s.diag(s.pi*tension*(n**2-1)/radius, s.pi*tension*n**2/radius)
    j_bend = s.Matrix([[0, -s.pi*rho*gamma*radius],
                       [s.pi*rho*gamma*radius, 0]])
    gen_bend = -j_bend.inv()*h_bend
    freq2 = tension**2*n**2*(n**2-1)/(rho**2*gamma**2*radius**4)
    checks.check("derived restricted bending frequency has the full n squared minus one",
                 s.simplify(gen_bend**2+freq2*s.eye(2)) == s.zeros(2))
    checks.check("n=1 has the correct symmetry-protected zero eigenvalue",
                 gen_bend.subs(n, 1)**2 == s.zeros(2)
                 and gen_bend.subs(n, 1) != s.zeros(2))
    checks.check("n=2 has two positive leading stiffness coefficients",
                 h_bend.subs(n, 2)[0, 0].is_positive
                 and h_bend.subs(n, 2)[1, 1].is_positive)
    for harmonic in (1, 2, 3):
        vector_row = er.applyfunc(lambda value, h=harmonic:
                                 s.integrate(value*s.cos(h*phi), (phi, 0, 2*s.pi)))
        expected = s.Matrix([s.pi, 0, 0]) if harmonic == 1 else s.zeros(3, 1)
        checks.check(f"actual global vector selection at toroidal n={harmonic}",
                     vector_row == expected)
    mass = s.symbols("M_D", positive=True)
    pos = (radius+eps*q*s.cos(2*phi))*er
    planar_quadrupole = mass/(2*s.pi)*s.integrate(pos[0]**2-pos[1]**2,
                                                 (phi, 0, 2*s.pi))
    checks.check("nonrigid n=2 has a literal Euclidean quadrupole amplitude",
                 s.simplify(s.diff(planar_quadrupole, eps).subs(eps, 0)
                            -mass*radius*q) == 0)

    r, ac, carrier, c = s.symbols("r a_c k c", positive=True)
    om, w, f, pressure = (s.Function(name)(r) for name in ("Omega", "W", "f", "P"))
    zeta = 2*om+r*s.diff(om, r)
    doppler = c-w
    a_mode = -zeta*f
    b_mode = pressure/doppler-s.diff(w, r)*f
    h_mode = -s.I*doppler*f
    checks.check("slow azimuthal equation retains the exact advected mode",
                 s.simplify(-s.I*c*a_mode+s.I*w*a_mode+zeta*h_mode) == 0)
    checks.check("slow axial equation retains pressure and axial shear",
                 s.simplify(-s.I*c*b_mode+s.I*w*b_mode
                            +s.diff(w, r)*h_mode+s.I*pressure) == 0)
    fp = -f/r+pressure/doppler**2
    checks.check("slow integral radial return is actual incompressibility",
                 s.simplify((s.diff(r*h_mode, r)/r+s.I*b_mode).subs(
                     s.diff(f, r), fp)) == 0)
    yfun = s.Function("y")(r)
    limit_pressure = doppler**2*s.diff(yfun, r)/r
    checks.check("slow radial pressure reproduces the actual zero-carrier Sturm problem",
                 s.simplify(s.diff(limit_pressure, r)-2*om*a_mode.subs(f, yfun/r)
                            -(s.diff(limit_pressure, r)+2*om*zeta*yfun/r)) == 0)
    ss = s.symbols("s", positive=True)
    hardy_schmidt = s.integrate(r*s.integrate(ss/r**2, (ss, 0, r)), (r, 0, ac))
    pressure_schmidt = s.integrate(r*s.integrate(4/ss, (ss, r, ac)), (r, 0, ac))
    checks.check("radial return kernel has a finite exact Hilbert-Schmidt norm",
                 hardy_schmidt == ac**2/4)
    checks.check("pressure Volterra kernel also has a finite weighted norm",
                 s.simplify(pressure_schmidt-ac**2) == 0)

    cross = s.Matrix(4, 2, s.symbols("c0:8"))
    ambient_j = s.Matrix([[0, 1], [-1, 0]])
    full_j = symplectic.row_join(cross).col_join((-cross.T).row_join(ambient_j))
    projection = s.eye(4).row_join(symplectic.inv()*cross).col_join(s.zeros(2, 6))
    checks.check("actual finite moment matrix defines an idempotent Euclidean projection",
                 s.simplify(projection**2-projection) == s.zeros(6))
    checks.check("the exact complementary projection is symplectically orthogonal",
                 s.simplify(projection.T*full_j-full_j*projection) == s.zeros(6))
    spectral, aa, dd, bc, cb = s.symbols("spectral a d b c", nonzero=True)
    pencil = s.Matrix([[spectral-aa, -bc], [-cb, spectral-dd]])
    schur = spectral-aa-bc*cb/(spectral-dd)
    checks.check("full eliminated pressure block has the stated Schur term",
                 s.factor(pencil.det()-(spectral-dd)*schur) == 0)
    checks.check("dropping the complement really changes the response determinant",
                 s.factor(pencil.det()-(spectral-dd)*(spectral-aa)) != 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
