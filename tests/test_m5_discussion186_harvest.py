"""Independent exact oracles for the bounded P255 discussion harvest."""

import itertools

import sympy as sp


def test_r17_pseudoscalar_reduces_to_three_r0_contractions():
    eta = (-1, 1, 1, 1)
    metric = sp.diag(*eta)
    jets = []
    for mu in range(4):
        a = sp.zeros(4)
        for i in range(4):
            for j in range(i, 4):
                a[i, j] = a[j, i] = sp.Symbol(f"a{mu}_{i}{j}")
        jets.append(a)
    curvature = {(m, n): jets[m] * metric * jets[n] - jets[n] * metric * jets[m]
                 for m in range(4) for n in range(4)}
    indices = itertools.product(range(4), repeat=4)
    i1 = sum(eta[m] * eta[n] * sp.trace(metric * curvature[m, n] * metric * curvature[m, n].T)
             for m in range(4) for n in range(m+1, 4))
    i2 = sum(curvature[m, n][a, b] * curvature[a, b][m, n] for m, n, a, b in indices)
    i3 = sum(eta[m] * eta[b] * curvature[m, n][a, b] * curvature[m, a][n, b]
             for m, n, a, b in itertools.product(range(4), repeat=4))
    x = sum(sp.LeviCivita(m, n, a, b) * eta[m] * eta[n] * curvature[m, n][a, b]
            for m, n, a, b in itertools.product(range(4), repeat=4)) / 2
    assert sp.expand(x*x + 2*i1 + i2 - 4*i3) == 0
    assert sp.expand(x*x + 2*i1 + i2 - 3*i3) != 0


def test_r13_compact_one_angle_has_flat_static_curvature_and_positive_clock():
    delta, psi = sp.symbols("delta psi", real=True)
    eta = sp.diag(-1, 1, 1, 1)
    rotation = sp.eye(4)
    rotation[1, 1] = rotation[2, 2] = sp.cos(psi)
    rotation[1, 2] = -sp.sin(psi)
    rotation[2, 1] = sp.sin(psi)
    m = rotation * sp.diag(8, 1, delta, 0) * rotation.T
    tangent = m.diff(psi)
    a_x, a_y = 2*tangent, 3*tangent
    assert (a_x*eta*a_y-a_y*eta*a_x).applyfunc(sp.simplify) == sp.zeros(4)
    other = sp.diag(0, 1, 0, 0)
    assert (a_x*eta*other-other*eta*a_x).subs({delta: sp.Rational(3, 10), psi: 0}) != sp.zeros(4)
    assert all(sp.trigsimp(sp.trace((m*eta)**p) - ((-8)**p+1+delta**p)) == 0
               for p in range(1, 5))
    generator = sp.zeros(4)
    generator[2, 3], generator[3, 2] = -1, 1
    clock = generator * m - m * generator
    f = clock * eta * tangent - tangent * eta * clock
    kinetic = sp.trigsimp(4 * sp.trace(eta * f * eta * f.T))
    expected = 8*(delta-1)**2*(1-(1-delta**2)*sp.cos(psi)**2)
    assert sp.trigsimp(kinetic-expected) == 0
    assert sp.simplify(kinetic.subs({delta: sp.Rational(3, 10), psi: 0})).is_positive
    assert sp.trigsimp(kinetic-2*expected) != 0


def test_r13_global_rotation_symmetry_does_not_preserve_fixed_exterior():
    eta = sp.diag(-1, 1, 1, 1)
    exterior = sp.diag(8, 1, sp.Rational(3, 10), 0)
    q = sp.eye(4)
    q[2, 2] = q[3, 3] = 0
    q[2, 3], q[3, 2] = -1, 1
    assert q*eta*q.T == eta
    assert q*exterior*q.T != exterior
    a = sp.diag(0, 1, 0, 0)
    b = sp.zeros(4)
    b[1, 2] = b[2, 1] = 1
    f = a*eta*b-b*eta*a
    fq = (q*a*q.T)*eta*(q*b*q.T)-(q*b*q.T)*eta*(q*a*q.T)
    assert f != sp.zeros(4)
    assert fq == q*f*q.T
    assert sp.trace(eta*f*eta*f.T) == sp.trace(eta*fq*eta*fq.T)


def test_r22_frame_gap_weights_are_adjacent_not_opposite():
    eigen = sp.symbols("lambda0:3", real=True)
    w = sp.symbols("w0:6", real=True)

    def hat(v):
        x, y, z = v
        return sp.Matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])

    d = sp.diag(*eigen)
    tangent = [hat(v)*d-d*hat(v) for v in (w[:3], w[3:])]
    comm = tangent[0]*tangent[1]-tangent[1]*tangent[0]
    cross = sp.Matrix(w[:3]).cross(sp.Matrix(w[3:]))
    gaps = [sp.prod(eigen[a]-eigen[b] for b in range(3) if b != a) for a in range(3)]
    assert sp.expand(4*sp.trace(comm*comm.T)-8*sum(gaps[a]**2*cross[a]**2 for a in range(3))) == 0
    opposite = sum((eigen[(a+1)%3]-eigen[(a+2)%3])**4*cross[a]**2 for a in range(3))
    assert sp.expand(4*sp.trace(comm*comm.T)-8*opposite) != 0
    assert [sp.sympify(g).subs(dict(zip(eigen, (1, sp.Rational(3, 10), sp.Rational(3, 10)))))
            for g in gaps] == [sp.Rational(49, 100), 0, 0]


def test_symmetric_four_eigenvalue_valley_slope_has_explicit_exclusions():
    nu, split, radius = sp.symbols("nu split radius", positive=True)
    values = (nu+split, nu, nu, nu)
    potential = (sum(x*x for x in values)-4*radius**2)**2
    derivative = sp.diff(potential, nu)
    center = {nu: radius, split: 0}
    assert derivative.subs(center) == 0
    assert -sp.diff(derivative, split).subs(center)/sp.diff(derivative, nu).subs(center) == -sp.Rational(1, 4)
    asymmetric = sum(c*(x-radius)**2 for c, x in zip((2, 1, 1, 1), values))
    assert -sp.diff(asymmetric, nu, split)/sp.diff(asymmetric, nu, 2) == -sp.Rational(2, 5)


def test_single_director_two_form_has_zero_wedge_but_two_forms_need_not():
    jets = sp.symbols("x0:8")

    def pullback(i, j):
        return jets[2*i]*jets[2*j+1]-jets[2*i+1]*jets[2*j]

    wedge = pullback(0, 1)*pullback(2, 3)-pullback(0, 2)*pullback(1, 3)+pullback(0, 3)*pullback(1, 2)
    assert sp.expand(wedge) == 0
    # H=u∧v and H'=s∧t from independent target fields. Each is decomposable;
    # their sum need not be: H01=1, H'23=1 gives (H+H')∧(H+H') != 0.
    first = {(0, 1): 1}
    second = {(2, 3): 1}
    both = {**first, **second}
    assert both.get((0, 1), 0)*both.get((2, 3), 0) - both.get((0, 2), 0)*both.get((1, 3), 0) + both.get((0, 3), 0)*both.get((1, 2), 0) == 1


def test_cp2_projector_curvature_needs_projected_gram_determinant():
    z = sp.Matrix([1, 0, 0])
    ui = sp.Matrix([0, 1+sp.I, 2-sp.I])
    uj = sp.Matrix([0, 2+3*sp.I, -1+sp.I])
    pi = ui*z.conjugate().T+z*ui.conjugate().T
    pj = uj*z.conjugate().T+z*uj.conjugate().T
    comm = pi*pj-pj*pi
    lhs = sp.simplify(sp.trace(comm*comm.conjugate().T))
    gij = sp.simplify((ui.conjugate().T*uj)[0])
    gram_det = (ui.conjugate().T*ui)[0]*(uj.conjugate().T*uj)[0]-sp.re(gij)**2
    flux = 2*sp.im(gij)
    assert sp.simplify(lhs-2*gram_det-sp.Rational(3, 2)*flux**2) == 0
    assert sp.simplify(lhs-sp.Rational(3, 2)*flux**2) != 0
