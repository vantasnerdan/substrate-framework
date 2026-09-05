"""Finite rational unequal-wave preparation: physical D energy and current."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger
from verify_stationary import transverse_law_average


def main():
    checks = CheckLedger("P251-0188-Chebyshev")
    aa = sp.Rational(1,100)
    energy = 1+aa**2
    psi = ef.add(ef.trig(1),ef.scale(ef.trig(2),aa))
    alpha = ef.add(ef.trig(1),ef.scale(ef.trig(2),-aa))
    u = (psi,ef.scale(ef.trig(2,kind="sin"),aa),ef.scale(ef.trig(1,kind="sin"),-1))

    def clean(field):
        return ef.add({wave:sp.cancel(value) for wave,value in field.items()})

    def add(*vectors):
        return tuple(clean(ef.add(*(vec[i] for vec in vectors))) for i in range(3))

    def scale(vector, value):
        return tuple(ef.scale(component,value) for component in vector)

    def pairing(left,right):
        return sp.cancel(ef.mul(left,right).get(ef.ZERO,0))

    def hh(field):
        return ef.add({wave:sum(q*q for q in wave)*value for wave,value in field.items()})

    kap = sp.Matrix(sp.symbols("kx ky kz",real=True))
    disp = sp.Matrix(sp.symbols("Dx Dy Dz",real=True))

    def law(polynomial):
        return transverse_law_average(polynomial,kap,disp)

    ss = sp.Symbol("s")
    harmonics = (3,5,7)
    weights = {n:sp.Rational(n*n*(n*n+1),2*(n*n-1)**2) for n in harmonics}
    denominator = sum(sp.Rational(n*n)/weights[n] for n in harmonics)
    coeffs = {n:2*n/(weights[n]*denominator) for n in harmonics}
    checks.check("finite odd harmonic weights satisfy the actual desingularized kernel moment", sum(n*coeffs[n] for n in harmonics) == 2)
    endpoint_cost = sum(weights[n]*coeffs[n]**2 for n in harmonics)/5
    checks.check("three exact harmonics lie below the positive energy-matching margin", endpoint_cost < sp.Rational(13,1280))

    def evaluate_poly(poly):
        result = {}
        power = {ef.ZERO:sp.Integer(1)}
        expanded = sp.Poly(poly,ss)
        for order in range(expanded.degree()+1):
            result = ef.add(result,ef.scale(power,expanded.nth(order)))
            power = ef.mul(power,psi)
        return clean(result)

    basis1 = psi
    basis3 = evaluate_poly(sp.chebyshevt(3,ss))
    fixed = evaluate_poly(coeffs[5]*sp.chebyshevt(5,ss)+coeffs[7]*sp.chebyshevt(7,ss))
    matrix = sp.Matrix([[pairing(ef.trig(1),basis1),pairing(ef.trig(1),basis3)],
                        [pairing(ef.trig(2),basis1),pairing(ef.trig(2),basis3)]])
    solved = matrix.inv()*sp.Matrix([-pairing(ef.trig(1),ef.add(alpha,fixed)),
                                    -pairing(ef.trig(2),ef.add(alpha,fixed))])
    first_integral = ef.add(fixed,ef.scale(basis1,solved[0]),ef.scale(basis3,solved[1]))
    rhs = clean(ef.add(alpha,first_integral))
    checks.check("finite positive contrast has a nonsingular exact moment matrix", matrix.det() != 0)
    checks.check("the actual finite polynomial forcing has zero kernel shell", all(sum(q*q for q in wave) != 1 for wave in rhs))
    phi = clean({wave:value/(sum(q*q for q in wave)-1) for wave,value in rhs.items()})
    z = (ef.add(ef.scale(alpha,-sp.Rational(1,2)),hh(phi)),ef.scale(ef.derivative(phi,2),-1),ef.derivative(phi,1))
    uy = (ef.trig(1),{},ef.scale(ef.trig(1,kind="sin"),-1))
    uz = (ef.scale(ef.trig(2),aa),ef.scale(ef.trig(2,kind="sin"),aa),{})
    force = ef.leray(ef.cross(uz,uy))
    generator = scale(ef.leray(add(ef.transport(u,z),ef.transport(z,u))),-1)
    checks.check("the degree-seven field solves the complete stationary Euler corrector", not any(add(generator,force)))
    d = kap[1]*disp[1]-kap[2]*disp[2]
    translation = tuple(ef.add(*(ef.scale(ef.derivative(u[i],j),-disp[j]) for j in range(3))) for i in range(3))
    kv = tuple({ef.ZERO:kap[i]} for i in range(3))
    lift = scale(ef.cross(kv,translation),-1)
    au = ef.add(*(ef.scale(u[i],kap[i]) for i in range(3)))
    ud = ef.add(*(ef.scale(u[i],disp[i]) for i in range(3)))
    rate = add(lift,scale(z,d),tuple(ef.scale(au,-disp[i]) for i in range(3)))
    current = tuple(ef.add(ef.scale(au,disp[i]),ef.scale(ud,kap[i])) for i in range(3))

    def minus(vector):
        shell = tuple(ef.add({wave:value for wave,value in component.items() if sum(q*q for q in wave) == 1}) for component in vector)
        return scale(add(ef.leray(shell),scale(ef.curl(shell),-1)),sp.Rational(1,2))

    amplitude = sp.Symbol("t",real=True)
    returned = add(scale(minus(rate),-1),scale(minus(current),amplitude))
    checks.check("the full correlated preparation return is a stationary Euler helicity field", not any(add(ef.curl(returned),returned)))
    final_rate = add(rate,returned)
    h2 = sp.factor(law(ef.inner(final_rate,final_rate)-pairing(au,au)*disp.dot(disp)))
    acc = sp.factor(law(pairing(au,au)*disp.dot(disp)+ef.inner(current,final_rate)))
    range_z = tuple({wave:value for wave,value in component.items() if sum(q*q for q in wave) > 1} for component in z)
    cost = sp.cancel(ef.inner(range_z,range_z)/5)
    print("finite rational range cost R=",cost)
    checks.check("finite rational microgeometry retains a strict positive matching window", cost < 13*energy/1280)
    checks.check("direct full-field energy gives the universal correlated kernel coefficient", sp.factor(h2-energy*amplitude**2/15+47*energy/240-cost) == 0)
    checks.check("direct full-current mean gives the physical restoring coefficient", sp.factor(acc-energy*(8*amplitude+13)/120) == 0)
    radicand = 100-960*cost/energy
    selected_amplitude = -(4+sp.sqrt(radicand))/8
    actual_a = energy*(sp.sqrt(radicand)-9)/120
    checks.check("the selected algebraic preparation has positive actual restoring energy", radicand > sp.Rational(361,4) and radicand < 100)
    checks.check("physical energy equals the actual mean response, with no supplied inertia", sp.simplify(h2.subs(amplitude,selected_amplitude)-actual_a) == 0 and sp.simplify(acc.subs(amplitude,selected_amplitude)+actual_a) == 0)
    checks.check("the actual micro return leaves the physical initial phase rho J unchanged", all(component.get(ef.ZERO,0) == 0 for component in returned))
    print("positive coefficient a=",actual_a)
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
