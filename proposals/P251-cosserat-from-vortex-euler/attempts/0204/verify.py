"""Exact geometry, packet, full-pressure jet and physical-density anchors."""

import sympy as s

from substrate_framework.euler_observation import (
    material_tag_fourier_dipole,
    material_tag_moments,
)
from substrate_framework.euler_phase import physical_scalar_chart
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0204-actual-ring-array-transfer")
    x, y, z, radius = s.symbols("x y z R", real=True, positive=True)
    coordinates = s.Matrix([x, y, z])
    tube = s.Matrix([(radius+x)*s.cos(z/radius)-radius,
                     y, (radius+x)*s.sin(z/radius)])
    jacobian = tube.jacobian(coordinates)
    determinant = s.simplify(jacobian.det())
    checks.check("actual tubular volume factor is retained",
                 s.simplify(determinant-1-x/radius) == 0)
    xi = s.Matrix([0, y-z, 2*x-z])
    pushed = jacobian*xi/determinant
    physical_div = s.trace(pushed.jacobian(coordinates)*jacobian.inv())
    checks.check("exact Piola preparation stays physically divergence-free",
                 s.simplify(s.trigsimp(physical_div)) == 0)
    checks.check("straightened core frame is not assigned unit volume at finite radius",
                 s.simplify(jacobian.T*jacobian
                            -s.diag(1, 1, (1+x/radius)**2)) == s.zeros(3))

    axial, carrier, length = s.symbols("q p L", real=True, positive=True)
    gaussian = length*s.exp(-length**2*(axial-carrier)**2/2)
    checks.check("finite packet amplitude has its actual Fourier normalization",
                 s.simplify(s.integrate(gaussian, (axial, -s.oo, s.oo))
                            /s.sqrt(2*s.pi)-1) == 0)
    checks.check("full packet phase mass scales by the squared spectral envelope",
                 s.simplify(s.integrate(gaussian**2, (axial, -s.oo, s.oo))
                            -s.sqrt(s.pi)*length) == 0)
    displacement, d3 = s.symbols("dp d3", real=True)
    cubic_average = d3*(displacement**3+3*displacement/length**2)
    checks.check("fiber current-two matching is not falsely exact after finite bandwidth",
                 s.diff(cubic_average, displacement).subs(displacement, 0)
                 == 3*d3/length**2)

    kk = s.symbols("k", real=True)
    kappa = s.Matrix(s.symbols("n0:3", real=True))
    phase = s.exp(-s.I*kk*kappa.dot(coordinates))
    potential = s.Matrix([y*z, x**2, x*z])

    def curl_local(vector):
        return s.Matrix([s.diff(vector[2], y)-s.diff(vector[1], z),
                         s.diff(vector[0], z)-s.diff(vector[2], x),
                         s.diff(vector[1], x)-s.diff(vector[0], y)])

    bloch_curl = curl_local(phase*potential)+s.I*kk*kappa.cross(phase*potential)
    checks.check("within-cell compensation is an actual Bloch-curl identity",
                 s.simplify(bloch_curl-phase*curl_local(potential)) == s.zeros(3, 1))
    qa, qb = s.symbols("xa xb", real=True)
    cross_phase = s.exp(-s.I*kk*qa)*s.exp(s.I*kk*qb)
    checks.check("cross-action spatial derivatives use separation not absolute radius",
                 s.simplify(s.diff(cross_phase, kk, 2).subs(kk, 0)+(qa-qb)**2) == 0)
    projector = s.eye(3)-kappa*kappa.T/kappa.dot(kappa)
    checks.check("actual zero-harmonic transverse pressure projector is idempotent",
                 s.simplify(projector**2-projector) == s.zeros(3))
    checks.check("actual zero-harmonic pressure projector removes only the longitudinal row",
                 s.simplify(projector*kappa) == s.zeros(3, 1))
    ntest = s.Matrix([1, 1, 0])/s.sqrt(2)
    ptest = s.eye(3)-ntest*ntest.T
    mean_force = -s.I*kk*ptest*s.diag(1, 2, 3)*ntest
    checks.check("retained anisotropic stress can create a nonzero Bloch mean jet",
                 s.diff(mean_force, kk) != s.zeros(3, 1))
    high_mean = projector*(kk**4*s.Matrix([1, 2, 3]))
    checks.check("high-moment preparation removes the initial singular mean through degree two",
                 all(s.diff(entry, kk, j).subs(kk, 0) == 0
                     for entry in high_mean for j in range(3)))
    laurent = s.symbols("w")
    checks.check("high angular harmonics kill full low Cartesian moments",
                 all(s.expand(laurent**12*((laurent+1/laurent)/2)**j
                              *(laurent+1/laurent)).coeff(laurent, 0) == 0
                     for j in range(5)))
    checks.check("low-harmonic mutation restores a full Cartesian moment",
                 s.expand(laurent*(laurent+1/laurent)).coeff(laurent, 0) == 1)

    time = s.symbols("t", real=True)
    coefficient, delta = s.symbols("J Delta", positive=True)
    f = s.Function("F")(time)
    g = s.Function("G")(time)
    row = s.Matrix([[f, g]])
    chart = physical_scalar_chart(
        coefficient*s.Matrix([[0, 1], [-1, 0]]), s.zeros(2), row,
        angle_rate=row.diff(time), angle_acceleration=row.diff(time, 2),
        generator_rate=s.zeros(2), spin=delta*row.diff(time),
    )
    wronskian = f*s.diff(g, time)-s.diff(f, time)*g
    numerator = s.diff(f, time)*s.diff(g, time, 2)-s.diff(f, time, 2)*s.diff(g, time)
    checks.check("actual moving observation mass uses the complete Wronskian",
                 s.simplify(chart.mass-coefficient/wronskian) == 0)
    checks.check("actual moving stiffness keeps both observation accelerations",
                 s.simplify(chart.stiffness-coefficient*numerator/wronskian**2) == 0)
    checks.check("measured quiet-tag spin is not renamed canonical inertia",
                 s.simplify(chart.spin_inertia-delta) == 0
                 and s.simplify(chart.spin_connection) == 0)

    wave = s.Matrix(s.symbols("K0:3", real=True))
    tag = material_tag_moments([1, 1], [[1, 0, 0], [-1, 0, 0]],
                               [[0, 1, 0], [0, -1, 0]])
    ambient = material_tag_moments([1, 1], [[2, 0, 0], [-2, 0, 0]],
                                   [[0, -s.Rational(1, 2), 0],
                                    [0, s.Rational(1, 2), 0]])
    tag_dipole = material_tag_fourier_dipole(wave, tag.spin, tag.shape_rate)
    ambient_dipole = material_tag_fourier_dipole(wave, ambient.spin, ambient.shape_rate)
    checks.check("complete point-current first moments can cancel with actual ambient return",
                 s.simplify(tag_dipole+ambient_dipole) == s.zeros(3, 1))
    checks.check("actual centroid-collapse current retains the nonzero tag dipole",
                 ambient_dipole != s.zeros(3, 1)
                 and s.simplify(ambient_dipole+tag_dipole) == s.zeros(3, 1))
    checks.check("the same exposing tag retains a nonzero symmetric shape-rate row",
                 tag.shape_rate != s.zeros(3))
    mcount = s.symbols("M", integer=True, positive=True)
    spacing, arcspacing, cube = s.symbols("d D c", positive=True)
    index = s.symbols("j", integer=True)
    count = mcount*s.summation(2*s.pi*(mcount*spacing+index*spacing)/arcspacing,
                               (index, 0, mcount-1))
    volume = (cube*mcount*spacing)**3
    checks.check("coaxial distributed-packet count has positive limiting density",
                 s.simplify(s.limit(count/volume, mcount, s.oo)
                            -3*s.pi/(cube**3*spacing**2*arcspacing)) == 0)
    checks.check("one densely prepared ring still dilutes in a radius-sized cube",
                 s.limit(radius/(cube*radius)**3, radius, s.oo) == 0)
    error, period, phase_mass = s.symbols("err P Mraw", positive=True)
    checks.check("relative physical action error uses one common volume normalization",
                 s.simplify((error/period**3)/(phase_mass/period**3)-error/phase_mass) == 0)
    checks.check("a fixed high band beats the explicit second-jet projector loss",
                 s.limit(period**(4-12), period, s.oo) == 0)
    checks.check("whole-law reconstruction does not multiply the raw phase density",
                 s.simplify((phase_mass/(3*period**3))/(delta/(3*period**3))
                            -phase_mass/delta) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
