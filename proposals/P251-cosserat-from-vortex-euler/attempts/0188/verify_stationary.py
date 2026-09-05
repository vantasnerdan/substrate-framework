"""Exact unequal-wave Euler corrector, mean/phase and action coefficient algebra."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def transverse_law_average(polynomial, kap, disp):
    """Exact Haar orthonormal-pair contraction of a quadratic D row."""
    conditional = 0
    for powers, coefficient in sp.Poly(sp.expand(polynomial), *disp).terms():
        indices = [i for i in range(3) for _ in range(powers[i])]
        if len(indices) != 2:
            raise ValueError("physical coefficient must be quadratic in polarization")
        i, j = indices
        conditional += coefficient * (sp.KroneckerDelta(i, j)-kap[i]*kap[j])/2
    result = 0
    for powers, coefficient in sp.Poly(sp.expand(conditional), *kap).terms():
        if any(power % 2 for power in powers):
            continue
        result += coefficient * sp.prod(sp.factorial2(power-1) for power in powers) / sp.factorial2(sum(powers)+1)
    return sp.simplify(result)


def main():
    checks = CheckLedger("P251-0188-stationary")
    a,b = sp.symbols("A B",positive=True)
    contrast = b*b-a*a
    psi = ef.add(ef.scale(ef.trig(1),b),ef.scale(ef.trig(2),a))
    alpha = ef.add(ef.scale(ef.trig(1),b),ef.scale(ef.trig(2),-a))
    u = (psi,ef.scale(ef.trig(2,kind="sin"),a),ef.scale(ef.trig(1,kind="sin"),-b))
    pressure = ef.scale(ef.add(*(ef.mul(component,component) for component in u)),-sp.Rational(1,2))
    gradp = tuple(ef.derivative(pressure,i) for i in range(3))

    def clean(field):
        return ef.add({wave:sp.cancel(value) for wave,value in field.items()})

    def add(*vectors):
        return tuple(clean(ef.add(*(vec[i] for vec in vectors))) for i in range(3))

    def scale(vector,value):
        return tuple(ef.scale(component,value) for component in vector)

    def hh(field):
        return ef.add({wave:sum(q*q for q in wave)*value for wave,value in field.items()})

    def transport(field):
        return ef.transport(u,(field,{},{}))[0]

    def pairing(left,right):
        return sp.cancel(ef.mul(left,right).get(ef.ZERO,0))

    checks.check("unequal field is exact stationary Euler with full pressure", not any(add(ef.transport(u,u),gradp)))
    checks.check("unequal field retains the actual curl eigenvalue", not any(add(ef.curl(u),u)))
    cubic = ef.mul(ef.mul(psi,psi),psi)
    moment_matrix = sp.Matrix([[pairing(ef.trig(1),psi),pairing(ef.trig(1),cubic)],
                               [pairing(ef.trig(2),psi),pairing(ef.trig(2),cubic)]])
    checks.check("first-integral moment determinant exposes the unequal-field repair", sp.factor(moment_matrix.det()-3*a*b*contrast/16) == 0)
    solved = sp.simplify(moment_matrix.inv()*sp.Matrix([-pairing(ef.trig(1),alpha),-pairing(ef.trig(2),alpha)]))
    rhs = clean(ef.add(alpha,ef.scale(psi,solved[0]),ef.scale(cubic,solved[1])))
    checks.check("solved actual forcing has zero first Fourier shell", all(sum(q*q for q in wave) != 1 for wave in rhs))
    phi = clean({wave:value/(sum(q*q for q in wave)-1) for wave,value in rhs.items()})
    expected = ef.add(ef.scale(ef.trig(1,3),b**3/(12*contrast)),ef.scale(ef.trig(2,3),a**3/(12*contrast)),
                      ef.scale(ef.mul(ef.trig(1,2),ef.trig(2)),b*b*a/contrast),
                      ef.scale(ef.mul(ef.trig(1),ef.trig(2,2)),b*a*a/contrast))
    checks.check("full scalar inversion gives the explicit finite Fourier corrector", not clean(ef.add(phi,ef.scale(expected,-1))))
    checks.check("stationary planar transport equation is exact", not clean(ef.add(transport(ef.add(hh(phi),ef.scale(phi,-1))),ef.scale(transport(alpha),-1))))
    r = ef.scale(alpha,-sp.Rational(1,2))
    correction = (ef.add(r,hh(phi)),ef.scale(ef.derivative(phi,2),-1),ef.derivative(phi,1))
    wave_y = (ef.scale(ef.trig(1),b),{},ef.scale(ef.trig(1,kind="sin"),-b))
    wave_z = (ef.scale(ef.trig(2),a),ef.scale(ef.trig(2,kind="sin"),a),{})
    force = ef.leray(ef.cross(wave_z,wave_y))
    velocity_generator = scale(ef.leray(add(ef.transport(u,correction),ef.transport(correction,u))),-1)
    checks.check("the full three-component stationary Euler corrector includes axial pressure", not any(add(velocity_generator,force)))
    checks.check("the constructed microscopic correction is mean-zero solenoidal", not ef.divergence(correction) and all(component.get(ef.ZERO,0) == 0 for component in correction))
    checks.check("equal amplitudes are an actual singular boundary, not a contradiction", moment_matrix.det().subs(b,a) == 0)

    kap = sp.Matrix(sp.symbols("kx ky kz",real=True))
    disp = sp.Matrix(sp.symbols("Dx Dy Dz",real=True))
    mx,my,mz = sp.symbols("mx my mz",real=True)
    mm = sp.diag(mx,my,mz)

    def law_average(polynomial):
        return transverse_law_average(polynomial,kap,disp)

    slow = sp.eye(3)-kap*kap.T
    ee = slow*mm*disp
    d = kap[1]*disp[1]-kap[2]*disp[2]
    de = kap[1]*ee[1]-kap[2]*ee[2]
    sigma = law_average(d*de)
    checks.check("actual projected tensor forcing includes its sixth moment", sp.simplify(sigma-(4*mx+19*(my+mz))/210) == 0)
    static = law_average((b*b*ee[1]*disp[1]+a*a*ee[2]*disp[2])/2)
    checks.check("the actual tensor static stress retains the transverse projection", sp.simplify(static-(b*b*(7*my+mx+my+mz)+a*a*(7*mz+mx+my+mz))/60) == 0)
    selected = {mx:sp.Rational(19,5),my:-sp.Rational(7,5),mz:sp.Rational(3,5),a:1,b:2}
    checks.check("a finite signed body tensor preserves mean trace and cancels forcing", (mx+my+mz).subs(selected) == 3 and sigma.subs(selected) == 0)
    checks.check("the projected unequal-field D comparison is genuinely restoring", static.subs(selected) == -sp.Rational(1,3))

    # Exact energy of the stationary-corrector class at one fixed rational
    # unequal field. This computes the coefficient; it does not select a sign.
    eta,xi,zeta = sp.symbols("eta xi zeta",real=True)
    ug = tuple({wave:value.subs({a:1,b:2}) for wave,value in component.items()} for component in u)
    psig = {wave:value.subs({a:1,b:2}) for wave,value in psi.items()}
    alphag = {wave:value.subs({a:1,b:2}) for wave,value in alpha.items()}
    phig = ef.add({wave:value.subs({a:1,b:2}) for wave,value in phi.items()},ef.scale(alphag,eta),ef.scale(psig,xi))
    rg = ef.add(ef.scale(alphag,-sp.Rational(1,2)),ef.scale(psig,zeta))
    zg = (ef.add(rg,hh(phig)),ef.scale(ef.derivative(phig,2),-1),ef.derivative(phig,1))
    translation = tuple(ef.add(*(ef.scale(ef.derivative(ug[i],j),-disp[j]) for j in range(3))) for i in range(3))
    kvector = tuple({ef.ZERO:kap[i]} for i in range(3))
    lift = scale(ef.cross(kvector,translation),-1)
    directional_u = ef.add(*(ef.scale(ug[i],kap[i]) for i in range(3)))
    rate = add(lift,scale(zg,d),tuple(ef.scale(directional_u,-disp[i]) for i in range(3)))
    # Enforce kappa.D=0 by the physical transverse contraction, not substitution
    # dividing by a selected component of kappa.
    h2 = sp.factor(law_average(ef.inner(rate,rate)-pairing(directional_u,directional_u)*disp.dot(disp)))
    restoring = (5*eta+3*xi)/10-sp.Rational(5,6)
    print("actual stationary D energy coefficient h2=",h2)
    print("actual restoring stress a=",restoring)
    print("full energy mismatch h2-a=",sp.factor(h2-restoring))
    checks.check("computed actual energy keeps quadratic preparation controls", sp.Poly(h2,eta,xi,zeta).total_degree() == 2)
    ee2,dd2 = sp.symbols("E Delta",positive=True)
    general_h2 = (ee2*(eta**2+xi**2)/5+ee2*zeta**2/10
                  +2*dd2*eta*xi/5+dd2*eta*zeta/5+ee2*xi*zeta/5
                  -ee2*eta/5-dd2*xi/5-dd2*zeta/5
                  +ee2*(25*ee2**2/(64*dd2**2)-sp.Rational(371,960)))
    checks.check("general contrast formula reproduces independently assembled full energy", sp.factor(general_h2.subs({ee2:5,dd2:3})-h2) == 0)
    mismatch = general_h2-(ee2*eta+dd2*xi)/10+ee2/6
    optimum = sp.solve([sp.diff(mismatch,var) for var in (eta,xi,zeta)],(eta,xi,zeta))
    minimum = sp.factor(mismatch.subs(optimum))
    checks.check("all stationary first-shell controls have an exact unique energy minimum", optimum == {eta:sp.Rational(3,4),xi:-dd2/(4*ee2),zeta:dd2/(2*ee2)})
    ratio = sp.symbols("r",positive=True)
    checks.check("every unequal contrast retains a strict energy mismatch", sp.factor(minimum.subs(dd2,ratio*ee2)-11*ee2/240-ee2*(1-ratio**2)*(375+12*ratio**2)/(960*ratio**2)) == 0)
    # A genuine velocity Jordan partner must satisfy first-integral moment
    # conditions, derived from H phi_D = -T B phi_return for a planar input.
    stationary_phi = clean(ef.add(phi,ef.scale(alpha,eta),ef.scale(psi,xi)))
    jordan_phi = ef.add(stationary_phi,ef.scale(alpha,-1))
    moment1 = pairing(psi,hh(jordan_phi))
    moment3 = pairing(cubic,hh(jordan_phi))
    jordan_controls = sp.solve([moment1,moment3],(eta,xi))
    print("necessary planar Jordan controls from psi and psi^3=",jordan_controls)
    fifth = ef.mul(cubic,ef.mul(psi,psi))
    fifth_residual = sp.factor(pairing(fifth,hh(jordan_phi)).subs(jordan_controls))
    print("remaining psi^5 Jordan moment=",fifth_residual)
    checks.check("Jordan moment controls are derived from the actual transport first integrals", sp.factor(moment1.subs(jordan_controls)) == 0 and sp.factor(moment3.subs(jordan_controls)) == 0)
    seventh = ef.mul(fifth,ef.mul(psi,psi))
    seventh_residual = sp.factor(pairing(seventh,hh(jordan_phi)).subs(jordan_controls))
    print("remaining psi^7 Jordan moment=",seventh_residual)
    numerator5 = sp.fraction(fifth_residual.subs(a,1))[0]
    numerator7 = sp.fraction(seventh_residual.subs(a,1))[0]
    common_factor = sp.gcd(numerator5,numerator7)
    checks.check("the polynomial stationary family has no full-field planar Jordan partner at any contrast", sp.Poly(common_factor,b).degree() == 0)
    tau = sp.symbols("tau",real=True)
    cubicg = {wave:value.subs({a:1,b:2}) for wave,value in cubic.items()}
    nonlinear_return = (ef.scale(cubicg,tau),{}, {})
    extended_rate = add(rate,scale(nonlinear_return,d))
    extended_h2 = sp.factor(law_average(ef.inner(extended_rate,extended_rate)-pairing(directional_u,directional_u)*disp.dot(disp)))
    nonlinear_velocity_generator = scale(ef.leray(add(ef.transport(ug,nonlinear_return),ef.transport(nonlinear_return,ug))),-1)
    checks.check("nonlinear axial first integral is an actual stationary Euler velocity return", not any(nonlinear_velocity_generator))
    cx = ef.add(ef.scale(directional_u,disp[0]),ef.scale(ef.add(*(ef.scale(ug[i],disp[i]) for i in range(3))),kap[0]))
    checks.check("the complete whole-law physical current of the nonlinear axial return vanishes", law_average(pairing(cx,ef.scale(cubicg,d))) == 0)
    print("actual nonlinear axial-return h2=",extended_h2)
    controls = (eta,xi,zeta,tau)
    nonlinear_minimizer = sp.solve([sp.diff(extended_h2-restoring,var) for var in controls],controls)
    print("nonlinear-return energy mismatch minimum=",sp.factor((extended_h2-restoring).subs(nonlinear_minimizer)))
    target = sp.symbols("a_target",real=True)
    # Solve the physical restoring row for eta before minimizing the remaining
    # true field controls; this keeps a positive response separate from a
    # possibly favorable unconstrained quadratic minimum.
    eta_for_target = sp.solve(restoring-target,eta)[0]
    constrained = sp.factor((extended_h2-restoring).subs(eta,eta_for_target))
    residual_minimizer = sp.solve([sp.diff(constrained,var) for var in (xi,zeta,tau)],(xi,zeta,tau))
    constrained_minimum = sp.factor(constrained.subs(residual_minimizer))
    print("nonlinear-return minimum mismatch at fixed actual restoring a=",constrained_minimum)
    checks.check("the extended positive-response test retains an exact quadratic mismatch", sp.Poly(constrained_minimum,target).degree() == 2)
    ud = ef.add(*(ef.scale(ug[i],disp[i]) for i in range(3)))
    full_c = tuple(ef.add(ef.scale(directional_u,disp[i]),ef.scale(ud,kap[i])) for i in range(3))
    projected_c = ef.leray(full_c)
    actual_r0 = sp.factor(law_average(pairing(directional_u,directional_u)*disp.dot(disp)+ef.inner(full_c,extended_rate)))
    print("DIRECT physical mean acceleration row at initial time=",actual_r0)
    print("actual projected current norm=",law_average(ef.inner(projected_c,projected_c)))
    checks.check("the declared restoring row agrees with the complete Euler/Lin initial observation", sp.factor(actual_r0+restoring) == 0)
    def minus_helicity(vector):
        shell = tuple(ef.add({wave:value for wave,value in component.items() if sum(q*q for q in wave) == 1}) for component in vector)
        return scale(add(ef.leray(shell),scale(ef.curl(shell),-1)),sp.Rational(1,2))

    zero_controls = {eta:0,xi:0,zeta:0,tau:0}
    base_rate = tuple(clean({wave:value.subs(zero_controls) for wave,value in component.items()}) for component in extended_rate)
    base_minus = minus_helicity(base_rate)
    current_minus = minus_helicity(projected_c)
    kernel_amplitude = sp.symbols("kernel_amplitude",real=True)
    kernel_return = add(scale(base_minus,-1),scale(current_minus,kernel_amplitude))
    checks.check("full correlated first-shell return has the actual stationary helicity", not any(add(ef.curl(kernel_return),kernel_return)))
    kernel_generator = scale(ef.leray(add(ef.transport(ug,kernel_return),ef.transport(kernel_return,ug))),-1)
    checks.check("correlated kernel cancellation solves the full Euler generator", not any(kernel_generator))
    kernel_rate = add(base_rate,kernel_return)
    kernel_h2 = sp.factor(law_average(ef.inner(kernel_rate,kernel_rate)-pairing(directional_u,directional_u)*disp.dot(disp)))
    kernel_acceleration = sp.factor(law_average(pairing(directional_u,directional_u)*disp.dot(disp)+ef.inner(full_c,kernel_rate)))
    print("full correlated-kernel energy=",kernel_h2)
    print("full correlated-kernel actual restoring row=",-kernel_acceleration)
    kernel_parameter = sp.solve(kernel_acceleration+target,kernel_amplitude)[0]
    kernel_mismatch = sp.factor((kernel_h2+kernel_acceleration).subs(kernel_amplitude,kernel_parameter))
    print("full correlated-kernel mismatch at fixed restoring a=",kernel_mismatch)
    range_velocity = tuple({wave:value for wave,value in component.items() if sum(q*q for q in wave) > 1} for component in base_rate)
    range_cost = law_average(ef.inner(range_velocity,range_velocity))
    print("actual cubic range cost=",range_cost)
    print("correlated-kernel mismatch without range energy=",sp.factor(kernel_mismatch-range_cost))
    checks.check("correlated kernel response is genuinely controllable", sp.diff(kernel_acceleration,kernel_amplitude) != 0)
    au_d = tuple(ef.scale(directional_u,disp[i]) for i in range(3))
    k_dot_translation = ef.add(*(ef.scale(translation[i],kap[i]) for i in range(3)))
    pressure_derivative = tuple(ef.add({wave:sp.I*wave[i]*value/sum(q*q for q in wave) for wave,value in k_dot_translation.items()}) for i in range(3))
    l1_constant = add(pressure_derivative,scale(ef.leray(au_d),-1))
    actual_b = add(l1_constant,scale(lift,-1))
    checks.check("velocity source uses the actual full Leray derivative and the same curl lift", not any(add(actual_b,ef.leray(add(au_d,lift)))))
    planar_input = {kap[0]:0,kap[1]:1,kap[2]:1,disp[0]:0,disp[1]:1,disp[2]:-1}
    planar_b = tuple(clean({wave:value.subs(planar_input) for wave,value in component.items()}) for component in actual_b)
    curl_b = ef.curl(planar_b)[0]
    stream_b = clean({wave:-value/sum(q*q for q in wave) for wave,value in curl_b.items()})
    checks.check("the exposing planar Jordan target includes its actual alpha source", not clean(ef.add(stream_b,ef.scale(alphag,-2))) and not clean(ef.add(planar_b[0],ef.scale(alphag,-1))))
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
