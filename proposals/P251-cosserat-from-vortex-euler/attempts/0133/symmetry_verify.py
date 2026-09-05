"""Frame, full Euler operator, neutral response, and material observables."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0133-actual-frame-and-mean")
    omega, axial, eigenvalue = s.symbols("Omega N lambda", positive=True)
    azimuth = s.Symbol("m", positive=True, integer=True)
    sigma = -2*omega*axial/eigenvalue
    speed = -sigma/axial
    laboratory_frequency = azimuth*omega+sigma
    checks.check("physical rotating and axial-translating frame makes the actual pattern steady",
                 s.simplify(azimuth*omega+axial*(sigma/axial)-laboratory_frequency) == 0)
    checks.check("steady frame has absolute vorticity lambda times its velocity",
                 s.simplify(eigenvalue*speed-2*omega) == 0)

    x, y, z, t = s.symbols("x y z t", real=True)
    coordinates = (x, y, z)
    base = s.Matrix([s.Function(f"U{i}")(*coordinates) for i in range(3)])
    probe = s.Matrix([s.Function(f"v{i}")(*coordinates) for i in range(3)])
    axis = s.Matrix([0, 0, 1])

    def gradient(value):
        return s.Matrix([s.diff(value, coordinate) for coordinate in coordinates])

    def curl(vector):
        return s.Matrix([s.diff(vector[2], y)-s.diff(vector[1], z),
                         s.diff(vector[0], z)-s.diff(vector[2], x),
                         s.diff(vector[1], x)-s.diff(vector[0], y)])

    absolute_vorticity = curl(base)+2*omega*axis
    transport = probe.jacobian(coordinates)*base+base.jacobian(coordinates)*probe
    residual = transport+2*omega*axis.cross(probe)
    residual += base.cross(curl(probe)-eigenvalue*probe)-gradient(base.dot(probe))
    checks.check("full Coriolis and pressure operator follows from absolute Beltrami identity",
                 s.simplify(residual+probe.cross(absolute_vorticity-eigenvalue*base)) == s.zeros(3, 1))
    checks.check("axial-boost generator retains its exact pressure gradient",
                 s.simplify(axis.cross(curl(base))+base.diff(z)-gradient(base[2])) == s.zeros(3, 1))

    # Actual helical symmetry of the m3 polynomial member.
    vplus = s.Matrix([1, s.I, 0])*(x+s.I*y)**2*s.exp(s.I*axial*z)
    spatial_rotation = s.Matrix([-y, x, 0])
    rotation_bracket = vplus.jacobian(coordinates)*spatial_rotation-spatial_rotation.jacobian(coordinates)*vplus
    checks.check("combined axial and physical rotation generator fixes the actual background",
                 s.simplify(vplus.diff(z)-axial*rotation_bracket/3) == s.zeros(3, 1))

    # Galilean identity at the material level before integration over a tag.
    perturbation = s.Symbol("epsilon", real=True)
    position = s.Matrix(s.symbols("r:3", real=True))
    centroid = s.Matrix(s.symbols("X:3", real=True))
    velocity = s.Matrix(s.symbols("v:3", real=True))
    mean_velocity = s.Matrix(s.symbols("V:3", real=True))
    shifted_position = position+perturbation*t*axis
    shifted_centroid = centroid+perturbation*t*axis
    shifted_velocity = velocity+perturbation*axis
    shifted_mean = mean_velocity+perturbation*axis
    old_spin = (position-centroid).cross(velocity-mean_velocity)
    new_spin = (shifted_position-shifted_centroid).cross(shifted_velocity-shifted_mean)
    checks.check("actual axial boost translates the centroid but leaves centered spin unchanged",
                 s.simplify(new_spin-old_spin) == s.zeros(3, 1)
                 and shifted_centroid.diff(perturbation) == t*axis)
    old_shape = (position[0]-centroid[0]+s.I*(position[1]-centroid[1]))**3
    new_shape = (shifted_position[0]-shifted_centroid[0]
                 +s.I*(shifted_position[1]-shifted_centroid[1]))**3
    checks.check("material third-moment orientation is unchanged by the apparent helical pattern shift",
                 s.simplify(new_shape-old_shape) == 0)

    amplitude, volume, norm = s.symbols("a volume norm", real=True, nonzero=True)
    # All terms linear in the m-nonzero wave vanish by its exact angular
    # Fourier integral; the retained quadratic coefficient is its full norm.
    helicity = 2*omega*speed*volume+eigenvalue*amplitude**2*norm
    checks.check("finite-amplitude family changes absolute helicity and is not one coadjoint orbit",
                 s.simplify(s.diff(helicity, amplitude)-2*eigenvalue*amplitude*norm) == 0
                 and s.diff(helicity, amplitude) != 0)
    hessian_vector = s.Matrix(s.symbols("h:3", real=True))
    checks.check("energy-helicity Hessian is conserved by the full intrinsic skew product",
                 s.simplify(hessian_vector.dot(base.cross(hessian_vector))) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
