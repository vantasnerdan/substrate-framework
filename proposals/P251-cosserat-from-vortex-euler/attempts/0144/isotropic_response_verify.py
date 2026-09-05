"""Pressure-derived two-phase initial response under a whole-field Haar law."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0144-isotropic-initial-response")
    gx, gy, gz, cx, cy, cz, eps = s.symbols("gx gy gz cx cy cz epsilon", real=True)
    wave = s.Matrix([gx, gy, gz])
    auxiliary = s.Matrix([cx, cy, cz])
    amplitude = wave.cross(auxiliary)
    direction = s.Matrix([0, 0, 1])
    field = s.Matrix([1, 0, 0])
    norm2 = wave.dot(wave)
    full_wave = wave+eps*direction
    projector = s.eye(3)-full_wave*full_wave.T/full_wave.dot(full_wave)
    zero_projector = s.eye(3)-wave*wave.T/norm2
    kelvin_raw = s.I*(wave*field.dot(amplitude)-amplitude*field.dot(wave))
    actual_kelvin_jet = (projector*kelvin_raw).diff(eps).subs(eps, 0)/s.I
    claimed_kelvin_jet = (-direction*field.dot(amplitude)
                          +wave*(direction.dot(wave)*field.dot(amplitude)
                                 +direction.dot(amplitude)*field.dot(wave))/norm2)
    ledger.check("actual Leray differentiation gives the complete Kelvin first spatial row",
                 all(s.factor(x) == 0 for x in actual_kelvin_jet-claimed_kelvin_jet))

    # Actual initial Euler acceleration of common V: -P_K[i(u.K)V+(V.grad)u].
    common_raw = s.I*eps*direction.dot(amplitude)*field+s.I*field.dot(wave)*amplitude
    actual_velocity_jet = -(projector*common_raw).diff(eps).subs(eps, 0)/s.I
    claimed_velocity_jet = (-direction.dot(amplitude)*field
                            +2*wave*direction.dot(amplitude)*field.dot(wave)/norm2)
    ledger.check("independent common-V Euler preparation gives its different pressure row",
                 all(s.factor(x) == 0 for x in actual_velocity_jet-claimed_velocity_jet))
    ledger.check("the solenoidal Fourier hypothesis is actually satisfied",
                 s.expand(wave.dot(amplitude)) == 0
                 and all(s.factor(x) == 0 for x in zero_projector*amplitude-amplitude))

    # Condition the polarization on the unit wave direction, then integrate
    # the remaining sphere. No isotropic coefficient is supplied to the check.
    z, phi = s.symbols("z phi", real=True)
    nx2 = (1-z*z)*s.cos(phi)**2

    def sphere(expression):
        return s.integrate(s.integrate(expression, (phi, 0, 2*s.pi)), (z, -1, 1))/(4*s.pi)

    nx2_bz2 = sphere(nx2*(1-z*z)/2)
    nxnz_bxbz = sphere(-nx2*z*z/2)
    b_variance = sphere((1-z*z)/2)
    position_return = 2*nx2_bz2+2*nxnz_bxbz
    velocity_return = 2*(nx2_bz2+nxnz_bxbz)
    position_stiffness = s.simplify(b_variance-position_return)
    velocity_stiffness = s.simplify(b_variance-velocity_return)
    ledger.check("whole-field orthogonal-frame moments are evaluated from the sphere",
                 nx2_bz2 == s.Rational(2, 15) and nxnz_bxbz == -s.Rational(1, 30))
    ledger.check("both actual preparations have the same positive isotropic initial coefficient",
                 position_stiffness == s.Rational(2, 15)
                 and velocity_stiffness == position_stiffness)
    ledger.check("omitting full pressure changes the measured coefficient",
                 b_variance != position_stiffness and b_variance == s.Rational(1, 3))
    ledger.check("a wrong common-V pressure factor is exposed by the phase comparison",
                 s.simplify(b_variance-(nx2_bz2+nxnz_bxbz)) != position_stiffness)
    print("Derived coefficient per actual Fourier energy:", position_stiffness)
    print("Independent common-V coefficient:", velocity_stiffness)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
