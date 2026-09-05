"""Exact algebra and independently differentiated smooth Euler fixtures."""

import pytest
import sympy as s

from substrate_framework import euler_fourier as ef


def background(harmonic=2):
    a, b = s.Rational(2, 3), s.Rational(3, 5)
    return (ef.scale(ef.trig(1, harmonic, "sin"), -b),
            ef.scale(ef.trig(0, harmonic, "sin"), a),
            ef.add(ef.scale(ef.trig(0, harmonic), a), ef.scale(ef.trig(1, harmonic), b)))


def generators(harmonic=2):
    return (({}, ef.trig(0, harmonic), ef.scale(ef.trig(0, harmonic, "sin"), -1)),
            (ef.trig(1, harmonic), {}, ef.trig(1, harmonic, "sin")))


def test_scalar_convolution_derivative_and_common_cell_average():
    cosine = ef.trig(0, s.Rational(3, 2))
    sine = ef.trig(0, s.Rational(3, 2), "sin")
    assert ef.derivative(cosine, 0) == ef.scale(sine, -s.Rational(3, 2))
    assert ef.add(ef.mul(cosine, cosine), ef.mul(sine, sine)) == {ef.ZERO: 1}
    assert ef.inner((cosine, {}, {}), (cosine, {}, {})) == s.Rational(1, 2)
    # Bilinear common-cell average is not a conjugating coefficient dot product.
    assert ef.inner(({(1, 0, 0): 1}, {}, {}), ({(1, 0, 0): 1}, {}, {})) == 0
    assert ef.trig(0, 0) == {ef.ZERO: 1}
    assert ef.trig(0, 0, "sin") == {}


def test_leray_mean_convention_idempotence_and_exact_gradient_removal():
    vector = ({ef.ZERO: s.Rational(2, 3), (1, 2, 0): 3}, {(1, 2, 0): 4}, {})
    projected = ef.leray(vector)
    assert ef.divergence(projected) == {}
    assert ef.leray(projected) == projected
    assert projected[0][ef.ZERO] == s.Rational(2, 3)
    assert ef.ZERO not in ef.leray(vector, mean_mode="zero")[0]
    scalar = ef.mul(ef.trig(0), ef.trig(1, 2, "sin"))
    assert ef.leray(tuple(ef.derivative(scalar, i) for i in range(3))) == ({}, {}, {})


@pytest.mark.parametrize("harmonic", [2, -2])
@pytest.mark.parametrize("rho", [s.Integer(1), s.Rational(5, 3)])
def test_complete_coadjoint_and_material_difference_with_density_and_lambda(harmonic, rho):
    u, xi = background(harmonic), generators(harmonic)
    velocities, H, omega = ef.coadjoint_matrices(u, xi, beltrami_eigenvalue=harmonic, density=rho)
    pressure = ef.scale(ef.add(*(ef.mul(component, component) for component in u)), -rho / 2)
    K = ef.material_jacobi_matrix(u, xi, pressure, density=rho)
    A = tuple(ef.material_kelvin_operator(u, field) for field in xi)
    gram = s.Matrix(2, 2, lambda i, j: rho * ef.inner(A[i], A[j]))
    assert H - K == gram
    assert gram != s.zeros(2)
    assert H == H.T and K == K.T and omega == -omega.T
    assert all(ef.divergence(v) == {} for v in velocities)
    # Actual Kelvin reconstruction is v-curl(xi cross u), not its negative.
    for i in range(2):
        curl_flux = ef.curl(ef.cross(xi[i], u))
        assert A[i] == tuple(ef.add(velocities[i][j], ef.scale(curl_flux[j], -1)) for j in range(3))
    # Independent signed constant mode: xi_0 cross omega averages lambda*a ex;
    # a curl has zero mean. A norm-only test would miss this reconstruction sign.
    assert A[0][0][ef.ZERO] == harmonic * s.Rational(2, 3)
    # A wrong helicity normalization or mean-zero projector changes this fixture.
    wrong_H = s.Matrix(2, 2, lambda i, j: rho * (ef.inner(velocities[i], velocities[j])
                       - ef.inner(velocities[i], ef.curl(velocities[j]))))
    assert wrong_H != H
    zero_mean_v = tuple(ef.leray(ef.cross(field, ef.curl(u)), mean_mode="zero") for field in xi)
    zero_mean_H = s.Matrix(2, 2, lambda i, j: rho * (ef.inner(zero_mean_v[i], zero_mean_v[j])
                           - ef.inner(zero_mean_v[i], ef.curl(zero_mean_v[j])) / harmonic))
    assert zero_mean_H - K != gram


def test_material_pressure_and_transport_against_independent_physical_calculus():
    x, y, z = s.symbols("x y z", real=True)
    coordinates = (x, y, z)
    rho = s.Rational(5, 3)
    u = s.Matrix([-s.Rational(3, 5) * s.sin(2 * y),
                  s.Rational(2, 3) * s.sin(2 * x),
                  s.Rational(2, 3) * s.cos(2 * x) + s.Rational(3, 5) * s.cos(2 * y)])
    xi = (s.Matrix([0, s.sin(2 * x), s.cos(2 * x)]),
          s.Matrix([s.sin(2 * y), 0, s.cos(2 * y)]))
    pressure = -rho * u.dot(u) / 2
    hessian = s.hessian(pressure, coordinates)
    convection = tuple(field.jacobian(coordinates) * u for field in xi)
    assert all(s.simplify(entry) == 0 for entry in u.jacobian(coordinates) * u
               + s.Matrix([s.diff(pressure, c) / rho for c in coordinates]))

    def mean(expression):
        # Independent real-space exact integrations, not Fourier zero-mode reuse.
        return s.simplify(s.integrate(s.integrate(s.expand(expression), (x, 0, 2 * s.pi)),
                                     (y, 0, 2 * s.pi)) / (2 * s.pi) ** 2)

    direct = s.Matrix(2, 2, lambda i, j: mean((xi[i].T * hessian * xi[j])[0]
                                             - rho * convection[i].dot(convection[j])))
    uf = background()
    pf = ef.scale(ef.add(*(ef.mul(c, c) for c in uf)), -rho / 2)
    xi_fourier = (({}, ef.trig(0, 2, "sin"), ef.trig(0, 2)),
                  (ef.trig(1, 2, "sin"), {}, ef.trig(1, 2)))
    actual = ef.material_jacobi_matrix(uf, xi_fourier, pf, density=rho)
    assert actual == direct
    assert actual != s.zeros(2)
    _, H, _ = ef.coadjoint_matrices(uf, xi_fourier, beltrami_eigenvalue=2, density=rho)
    A = tuple(ef.material_kelvin_operator(uf, field) for field in xi_fourier)
    assert H - direct == s.Matrix(2, 2, lambda i, j: rho * ef.inner(A[i], A[j]))
    missing_pressure_density = ef.material_jacobi_matrix(uf, xi_fourier, ef.scale(pf, 1 / rho), density=rho)
    assert missing_pressure_density != actual


def test_nonzero_kks_has_density_and_time_reversal_sign():
    u = background(1)
    axial = ef.trig(2)
    angle = (ef.scale(ef.mul(ef.trig(1, kind="sin"), axial), -1),
             ef.mul(ef.trig(0, kind="sin"), axial), {})
    shape = (ef.mul(ef.mul(ef.trig(0, kind="sin"), ef.trig(1)), axial),
             ef.scale(ef.mul(ef.mul(ef.trig(0), ef.trig(1, kind="sin")), axial), -1), {})
    _, H, omega = ef.coadjoint_matrices(u, (angle, shape), beltrami_eigenvalue=1, density=3)
    _, Hminus, minus = ef.coadjoint_matrices(tuple(ef.scale(c, -1) for c in u), (angle, shape), beltrami_eigenvalue=1, density=3)
    _, Hone, one = ef.coadjoint_matrices(u, (angle, shape), beltrami_eigenvalue=1)
    assert omega != s.zeros(2)
    assert Hminus == H and minus == -omega
    assert H == 3 * Hone and omega == 3 * one


def test_invalid_domains_are_rejected_without_silent_convention_changes():
    with pytest.raises(ValueError, match="nonzero"):
        ef.coadjoint_matrices(background(), generators(), beltrami_eigenvalue=0)
    with pytest.raises(ValueError, match="background"):
        ef.coadjoint_matrices(background(), generators(), beltrami_eigenvalue=1)
    with pytest.raises(ValueError, match="density"):
        ef.coadjoint_matrices(background(), generators(), beltrami_eigenvalue=2, density=-1)
    with pytest.raises(ValueError, match="generators"):
        ef.coadjoint_matrices(background(), ((ef.trig(0), {}, {}),), beltrami_eigenvalue=2)
    with pytest.raises(ValueError, match="mean_mode"):
        ef.leray(({}, {}, {}), mean_mode="invented")
    with pytest.raises(ValueError, match="three"):
        ef.curl(({}, {}))
    with pytest.raises(ValueError, match="kind"):
        ef.trig(0, kind="tan")
    with pytest.raises(ValueError, match="rational"):
        ef.trig(0, s.sqrt(2))


@pytest.mark.parametrize("density", [s.I, 1 + s.I, s.oo, -s.oo, s.zoo, s.nan])
def test_explicit_complex_or_nonfinite_density_is_rejected(density):
    with pytest.raises(ValueError, match="density"):
        ef.coadjoint_matrices(background(), generators(), beltrami_eigenvalue=2, density=density)
    with pytest.raises(ValueError, match="density"):
        ef.material_jacobi_matrix(background(), generators(), {}, density=density)


def test_undetermined_density_symbol_retains_explicit_positive_hypothesis():
    rho = s.Symbol("rho")
    _, H, _ = ef.coadjoint_matrices(background(), generators(), beltrami_eigenvalue=2, density=rho)
    _, Hunit, _ = ef.coadjoint_matrices(background(), generators(), beltrami_eigenvalue=2)
    assert H == rho * Hunit
