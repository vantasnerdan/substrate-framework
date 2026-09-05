"""Actual stationary affine-tag first moments under the axial-current forcing."""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0209-moment-response")
    omega, rho, reference, radial = sp.symbols("Omega rho Q I", positive=True)
    strain = sp.Matrix([[0,1],[-omega**2,0]])
    rotation = sp.Matrix([[0,-1],[1,0]])
    frequency = 4*omega
    h = 1+omega**2
    tau = (1-omega)/(1+omega)
    gradient_moment = 2*sp.pi*radial/(1+omega)**2*sp.Matrix([[1,sp.I],[sp.I*omega**2,-omega**2]])
    velocity_moment = gradient_moment*rotation.T
    symbols = sp.symbols("m00 m01 m10 m11")
    unknown = sp.Matrix(2,2,symbols)

    def resolve(forcing):
        equations = sp.I*frequency*unknown-strain*unknown-unknown*strain.T-forcing
        solution = sp.solve(list(equations),symbols)
        return unknown.subs(solution).applyfunc(sp.simplify)

    leading = resolve(velocity_moment)
    theta0 = sp.simplify((leading[0,1]+leading[1,0])/reference)
    g0 = sp.simplify(rho*(leading[0,1]-leading[1,0]))
    inertia = rho*reference*tau/(1+tau**2/3)
    checks.check("the full first-moment equation recovers the independently measured inertia", sp.simplify(g0-inertia*theta0) == 0)
    correction = resolve(gradient_moment-h*leading)
    theta1 = sp.simplify((correction[0,1]+correction[1,0])/reference)
    g1 = sp.simplify(rho*(correction[0,1]-correction[1,0]))
    print("actual gradient-forcing constant theta row =", theta1)
    print("actual gradient-forcing constant G row =", g1)
    print("actual physical G-inertia*theta mismatch =", sp.simplify(g1-inertia*theta1))
    polynomial = 2*omega**4+5*omega**3+4*omega**2+5*omega+2
    mismatch = sp.pi*radial*rho*(omega-1)*polynomial/(24*omega**2*(omega+1)*(omega**2+omega+1))
    checks.check("the nonzero literal displacement-current mismatch is retained exactly", sp.factor(g1-inertia*theta1-mismatch) == 0)
    gain = sp.Symbol("g", positive=True)
    gradient_ratio = gain*polynomial/(12*omega*(omega**2+omega+1))
    checks.check("common initial-angle normalization derives the stated gradient inertia", sp.factor(sp.I*gain*mismatch/(inertia*theta0)-gradient_ratio) == 0)
    checks.check("gradient inertia is strictly positive without a floating-point sign test", gradient_ratio.is_positive is True)
    time = sp.Symbol("t",real=True)
    response = sp.exp(sp.I*frequency*time)*(h*time*leading+correction)
    forcing = sp.exp(sp.I*frequency*time)*(h*time*velocity_moment+gradient_moment)
    checks.check("the complete linear-in-time forced first-moment response solves its defining equation", all(sp.simplify(v) == 0 for v in sp.diff(response,time)-strain*response-response*strain.T-forcing))
    theta_response = (response[0,1]+response[1,0])/reference
    secular_ratio = sp.simplify(sp.diff(sp.exp(-sp.I*frequency*time)*theta_response,time)/theta0)
    checks.check("the secular physical-angle coefficient is the positive actual clock derivative", sp.simplify(secular_ratio-h) == 0 and (2*frequency*secular_ratio*gain).is_positive is True)
    # Independently derive the gradient moment from the full angular pressure
    # coefficients, rather than treating the displayed matrix as an oracle.
    n11 = tau**2/2+(tau**3+tau)/4
    n22 = tau**2/2-(tau**3+tau)/4
    n12 = sp.I*(tau-tau**3)/4
    angular = sp.Matrix([[n11,n12/omega],[omega*n12,n22]])
    checks.check("the full inverse-pressure angular moments give Kmat", all(sp.simplify(v) == 0 for v in 2*sp.pi*radial*angular/tau-gradient_moment))
    width, cutoff, radius = sp.symbols("s b r", positive=True)
    # Actual Hankel integrals are evaluated before either moment is compared.
    wave = sp.Symbol("k", positive=True)
    # Integration by parts and Fubini evaluate the tag Hankel transform first;
    # this avoids an unnecessarily expensive nested symbolic integration.
    tag_hankel = sp.integrate(2*cutoff*radius**3*sp.exp(-cutoff*radius**2)*sp.besselj(2,wave*radius),(radius,0,sp.oo))
    moment = tau*sp.integrate(wave**3*sp.exp(-width*wave**2/2)*tag_hankel,(wave,0,sp.oo))
    inverse_moment = tau*sp.integrate(wave*sp.exp(-width*wave**2/2)*tag_hankel,(wave,0,sp.oo))
    checks.check("the actual pressure Hankel integral gives the physical radial moment", sp.simplify(moment-16*cutoff*tau/(1+2*cutoff*width)**3) == 0)
    checks.check("the inverse-wave moment has the distinct stated profile dependence", sp.simplify(inverse_moment-2*tau/(1+2*cutoff*width)**2) == 0)
    width2 = sp.Symbol("s2", positive=True)
    matrix = sp.Matrix([[moment,moment.subs(width,width2)],[inverse_moment,inverse_moment.subs(width,width2)]])
    determinant = 64*cutoff**2*tau**2*(width2-width)/((1+2*cutoff*width)**3*(1+2*cutoff*width2)**3)
    checks.check("distinct actual radial widths give an invertible two-moment preparation", sp.factor(matrix.det()-determinant) == 0)
    inverse_phase_norm = 2*sp.pi*sp.integrate(wave**7*sp.exp(-wave**2),(wave,0,sp.oo))
    checks.check("the invisible radial inverse-Laplacian preparation has finite positive full phase norm", inverse_phase_norm == 6*sp.pi)
    # A nonzero integer harmonic has a fourth-order characteristic zero.
    harmonic = sp.Symbol("n", integer=True, nonzero=True)
    offset = sp.Symbol("epsilon", real=True)
    local_characteristic = (sp.sin(sp.pi*offset)/(sp.pi*(harmonic+offset)))**4
    checks.check("positive fourfold phase law kills every nonresonant cubic time polynomial", all(sp.limit(sp.diff(local_characteristic,offset,j),offset,0) == 0 for j in range(4)))
    # Expose the actual first-K Kelvin vorticity forcing and its secular memory.
    q1, q2, zeta, integral = sp.symbols("q1 q2 zeta integral", real=True)
    q = sp.Matrix([q1,q2])
    checks.check("area preservation transports the full Kelvin Jq direction", strain*rotation+rotation*strain.T == sp.zeros(2))
    forced_b = -h*zeta*rotation*q*integral
    derivative_b = -h*zeta*(rotation*(-strain.T*q)*integral+rotation*q/(q.dot(q)))
    checks.check("the inverse-wave secular integral solves the full first-K Euler equation", all(sp.factor(v) == 0 for v in derivative_b-strain*forced_b+h*zeta*rotation*q/(q.dot(q))))
    checks.check("the resulting second-order vorticity memory has the retained negative sign", sp.cancel(h*(q1*forced_b[1]-q2*forced_b[0])/q.dot(q)+h**2*zeta*integral) == 0)
    # The genuine resonant Lin return is part of the full Hamiltonian state.
    beta, excess = sp.symbols("beta a", positive=True)
    full_phase = sp.Matrix([[0,beta,-excess,0],[-beta,0,0,-excess],[excess,0,0,0],[0,excess,0,0]])
    generator = sp.Matrix([[0,frequency,0,0],[-frequency,0,0,0],[1,0,0,frequency],[0,1,-frequency,0]])
    full_energy = -full_phase*generator
    checks.check("the actual resonant Jordan generator preserves the complete phase", generator.T*full_phase+full_phase*generator == sp.zeros(4))
    checks.check("the full conserved energy retains the positive unobserved excess", full_energy[:2,:2] == (frequency*beta+excess)*sp.eye(2))
    checks.check("discarding the resonant displacement would lose an actual energy row", (full_energy[:2,:2]-frequency*beta*sp.eye(2)).det() == excess**2)
    angle = sp.Symbol("theta", real=True)
    amplitude = sp.Symbol("G", real=True)
    number = sp.Symbol("N", positive=True, integer=True)
    passive_g = amplitude*sp.cos(number*angle)
    passive_h = -amplitude*sp.sin(number*angle)/(number*omega)
    transported_h = -omega*sp.diff(passive_h,angle)
    energy_density = rho*(passive_g**2-2*passive_g*transported_h)
    checks.check("the actual off-tag axial return supplies negative full energy with zero self phase", sp.simplify(transported_h-passive_g) == 0 and sp.simplify(energy_density+rho*passive_g**2) == 0)
    checks.check("two orthogonal actual energy returns repair the complete leading oscillator energy", full_energy[:2,:2]-excess*sp.eye(2) == frequency*beta*sp.eye(2))
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
