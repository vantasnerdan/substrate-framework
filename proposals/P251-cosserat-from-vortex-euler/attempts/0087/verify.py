"""Exact mass-preserving mean, spin-current and varied kinetic identities."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0087-volume-preserving-GLM")
    x, y, z, tau = s.symbols("x y z tau", real=True)
    coords = (x, y, z)
    amplitude = s.Symbol("d", real=True)
    pushed_density = 1+amplitude*s.cos(x)
    interpolation = (1-tau)*pushed_density+tau
    poisson = -amplitude*s.cos(x)
    correcting_velocity = s.diff(poisson, x)/interpolation
    ledger.check("physical-space Poisson correction has the actual density source",
                 s.simplify(s.diff(poisson, x, 2)-pushed_density+1) == 0)
    ledger.check("auxiliary mean flow pushes arithmetic-mean density to uniform volume exactly",
                 s.simplify(s.diff(interpolation, tau)
                            +s.diff(interpolation*correcting_velocity, x)) == 0)
    ledger.check("a volume-preserving deterministic input needs no correction",
                 s.simplify(correcting_velocity.subs(amplitude, 0)) == 0)

    def gradient(f):
        return s.Matrix([s.diff(f, c) for c in coords])

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    a, b, adot, bdot, epsilon = s.symbols("a b adot bdot epsilon", real=True)
    xi = s.Matrix([a*s.sin(y), b*s.sin(x), 0])
    xi_t = s.Matrix([adot*s.sin(y), bdot*s.sin(x), 0])
    jac = xi.jacobian(coords)
    advective = jac*xi
    covariance = xi*xi.T
    divergence_covariance = s.Matrix([sum(s.diff(covariance[i, j], coords[j])
                                         for j in range(3)) for i in range(3)])
    ledger.check("actual fluctuation generators are divergence free",
                 s.trace(jac) == 0)
    ledger.check("arithmetic material-mean correction comes from the actual covariance divergence",
                 s.simplify(advective-divergence_covariance) == s.zeros(3, 1))
    ledger.check("explicit fluctuation example has a purely longitudinal mean-volume correction",
                 s.simplify(advective-gradient(-a*b*s.cos(x)*s.cos(y))) == s.zeros(3, 1))
    second_map = s.Matrix(coords)+epsilon*xi+epsilon**2*advective/2
    determinant = s.expand(second_map.jacobian(coords).det())
    ledger.check("flow-map construction preserves volume through the claimed amplitude order",
                 s.simplify(determinant.coeff(epsilon, 1)) == 0
                 and s.simplify(determinant.coeff(epsilon, 2)) == 0)
    physical_mean2 = (jac*xi_t-xi_t.jacobian(coords)*xi)/2
    spin = xi.cross(xi_t)
    ledger.check("direct inverse-flow velocity gives half the fluctuation-spin curl",
                 s.simplify(physical_mean2-curl(spin)/2) == s.zeros(3, 1))
    ledger.check("spin-current shift is nonzero although the corrected material mean is zero",
                 s.simplify(physical_mean2) != s.zeros(3, 1))
    ledger.check("wrong spin-current sign fails the actual velocity expansion",
                 s.simplify(physical_mean2+curl(spin)/2) != s.zeros(3, 1))

    # Full same-material kinetic pullback, without setting pseudomomentum zero.
    deformation = s.Matrix([[1, a, 0], [0, 1, b], [0, 0, 1]])
    shape_columns = s.Matrix([[1, 0], [0, 1], [1, 1]])
    mean_velocity = s.Matrix(s.symbols("U0:3", real=True))
    rates = s.Matrix(s.symbols("qdot0:2", real=True))
    rho = s.Symbol("rho", positive=True)
    actual_velocity = deformation*mean_velocity+shape_columns*rates
    kinetic = rho*actual_velocity.dot(actual_velocity)/2
    mean_momentum = s.Matrix([s.diff(kinetic, component) for component in mean_velocity])
    internal_momentum = s.Matrix([s.diff(kinetic, component) for component in rates])
    ledger.check("full pullback gives the mean one-form momentum, not an assumed rho times velocity",
                 s.simplify(mean_momentum-rho*deformation.T*actual_velocity) == s.zeros(3, 1))
    ledger.check("retaining fluctuation variations gives the independent internal momentum equation",
                 s.simplify(internal_momentum-rho*shape_columns.T*actual_velocity) == s.zeros(2, 1))
    ledger.check("pseudomomentum mutation would remove nonzero same-action terms",
                 s.simplify(mean_momentum-rho*mean_velocity) != s.zeros(3, 1))

    # A true rotating material covariance gives the spin/metric coefficient.
    rate = s.Matrix(s.symbols("Omega0:3", real=True))
    c = s.Symbol("c", positive=True)
    covariance = c*s.eye(3)
    spin_coefficient = rho*(s.trace(covariance)*s.eye(3)-covariance)
    ledger.check("actual isotropic rotating-displacement spin equals its geometric kinetic metric",
                 spin_coefficient*rate == 2*rho*c*rate)
    # Projected rotation is divergence-free in a Fourier mode even when the
    # polarization projector and rigid-vector rotation do not commute.
    wave = s.Matrix([1, 1, 1])
    projector = s.eye(3)-wave*wave.T/3
    jz = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    eta0 = projector*s.Matrix(s.symbols("eta0:3", real=True))
    tangent = projector*jz*eta0
    ledger.check("projected material rotation gives spin=kinetic norm without assuming projection commutes",
                 s.simplify((jz*eta0).dot(tangent)-tangent.dot(tangent)) == 0
                 and projector*jz != jz*projector)
    ledger.check("projected rotational fluctuation evolution preserves its norm",
                 s.simplify(eta0.dot(tangent)) == 0)
    scalar, scalar_dot = s.symbols("q qdot", real=True)
    one_profile = s.Matrix(s.symbols("Y0:3", real=True))
    ledger.check("a positive one-profile oscillation alone does not create fluctuation spin",
                 (scalar*one_profile).cross(scalar_dot*one_profile) == s.zeros(3, 1))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
