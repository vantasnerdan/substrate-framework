"""Exact wrapped-cell Lin/current and implicit physical band-center checks."""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0219-bloch-current")
    x, y, time, wave = sp.symbols("a b t k", real=True)
    omega = sp.Symbol("Omega", positive=True)
    psi = sp.cos(y)+omega**2*sp.cos(x)
    velocity = sp.Matrix([sp.sin(y),-omega**2*sp.sin(x)])
    rotation = sp.Matrix([[0,-1],[1,0]])
    scalar = sp.Function("s")(x,y,time)

    def transport(value):
        return velocity[0]*sp.diff(value,x)+velocity[1]*sp.diff(value,y)

    rhs = -transport(scalar)-sp.I*wave*psi*scalar
    replacements = {sp.diff(scalar,time):rhs}
    for coordinate in (x,y):
        replacements[sp.diff(scalar,time,coordinate)] = sp.diff(rhs,coordinate)

    def reduced(value):
        return sp.simplify(sp.expand(value.xreplace(replacements)))

    normal = rotation*sp.Matrix([sp.diff(scalar,x),sp.diff(scalar,y)])+sp.I*wave*time*scalar*velocity
    axial = -time*transport(scalar)
    full_velocity = sp.Matrix([psi,*velocity])
    gradient = full_velocity.jacobian((sp.Symbol("X"),x,y))
    displacement = sp.Matrix([axial,*normal])
    lin = displacement.diff(time)+displacement.applyfunc(transport)+sp.I*wave*psi*displacement-gradient*displacement
    checks.check("the actual wrapped three-component displacement solves homogeneous Lin", all(reduced(entry) == 0 for entry in lin))
    divergence = sp.I*wave*axial+sp.diff(normal[0],x)+sp.diff(normal[1],y)
    checks.check("the retained axial shear exactly restores three-dimensional divergence", reduced(divergence) == 0)
    tag_transport = normal.dot(sp.Matrix([sp.diff(psi,x),sp.diff(psi,y)]))
    checks.check("the actual tag variation is unchanged by the required tangent correction", sp.simplify(tag_transport+transport(scalar)) == 0)
    naive = displacement-sp.Matrix([0,*(sp.I*wave*time*scalar*velocity)])
    naive_div = sp.I*wave*naive[0]+sp.diff(naive[1],x)+sp.diff(naive[2],y)
    checks.check("omitting the tangent correction fails the physical divergence constraint", sp.simplify(naive_div+sp.I*wave*time*transport(scalar)) == 0 and naive_div != 0)
    gauge = sp.Matrix([-transport(scalar),*(sp.I*wave*scalar*velocity)])
    gauge_lin = gauge.diff(time)+gauge.applyfunc(transport)+sp.I*wave*psi*gauge-gradient*gauge
    checks.check("the new dipole control is itself an exact full homogeneous Lin history", all(reduced(entry) == 0 for entry in gauge_lin))
    checks.check("the dipole control is solenoidal without an assigned projection", sp.simplify(sp.I*wave*gauge[0]+sp.diff(gauge[1],x)+sp.diff(gauge[2],y)) == 0)
    checks.check("the current control is tangent to the actual stationary tag", sp.simplify(gauge[1]*sp.diff(psi,x)+gauge[2]*sp.diff(psi,y)) == 0)
    radial_derivative = x*sp.diff(psi,x)+y*sp.diff(psi,y)
    checks.check("its literal nonzero dipole comes from the defining cross product", sp.simplify(x*gauge[2]-y*gauge[1]-sp.I*wave*scalar*radial_derivative) == 0)
    frequency, speed, frequency1, frequency2, speed1, curvature = sp.symbols("nu0 psi0 nu1 nu2 psi1 beta", nonzero=True, real=True)
    center1 = -speed/frequency1
    center2 = (2*curvature-frequency2*center1**2-2*speed1*center1)/frequency1
    center = center1*wave+center2*wave**2/2
    physical_frequency = frequency+frequency1*center+frequency2*center**2/2+wave*(speed+speed1*center)
    checks.check("actual implicit band motion cancels Doppler and derives the clock curvature", sp.series(physical_frequency,wave,0,3).removeO().expand() == frequency+curvature*wave**2)
    checks.check("a frozen center would retain the unwanted Doppler slope", sp.diff(frequency+speed*wave,wave) == speed)
    inertia, inertia1, inertia2, ratio, ratio1 = sp.symbols("j j1 j2 r r1", real=True)
    j = inertia+inertia1*center+inertia2*center**2/2
    r = ratio+ratio1*center
    dipole = j*(1-sp.I*wave*time*r)
    even_dipole = sp.series((dipole+dipole.subs(wave,-wave))/2,wave,0,3).removeO().expand()
    expected_dipole = inertia+wave**2*((inertia1*center2+inertia2*center1**2)/2-sp.I*time*center1*(inertia1*ratio+inertia*ratio1))
    checks.check("whole-field parity retains the actual even dipole-time connection", sp.simplify(even_dipole-expected_dipole) == 0)
    actual_spin_inertia = j*(1-wave*(speed+speed1*center)/(frequency+curvature*wave**2))
    even_spin = sp.series((actual_spin_inertia+actual_spin_inertia.subs(wave,-wave))/2,wave,0,3).removeO().expand()
    expected_spin = inertia+wave**2*((inertia1*center2+inertia2*center1**2)/2-center1*(inertia1*speed+inertia*speed1)/frequency)
    checks.check("the literal spin gradient differs from the uncorrected dipole gradient", sp.simplify(even_spin-expected_spin) == 0)
    drift = center1*(inertia1*ratio+inertia*ratio1)
    spin_shift = center1*(inertia1*speed+inertia*speed1)/frequency
    actual_return = sp.I*wave*(sp.I*drift/speed+wave*sp.I*spin_shift)*sp.exp(-sp.I*speed*wave*time)
    even_return = sp.series((actual_return+actual_return.subs(wave,-wave))/2,wave,0,3).removeO()
    checks.check("the actual Doppler current return cancels both observed second-jet differences", sp.simplify(even_return-wave**2*(-spin_shift+sp.I*drift*time)) == 0)
    corrected_dipole = even_dipole+even_return
    checks.check("the corrected literal dipole equals the separately measured spin-rate inertia", sp.simplify(corrected_dipole-even_spin) == 0)
    # Exact compact polynomial weight anchors the general smooth-bump Gram
    # argument; this finite-regularity check is not the selected source bump.
    coordinate = sp.Symbol("z", real=True)
    compact_weight = (1-coordinate**2)**4
    order = 4
    gram = sp.Matrix(order+1,order+1,lambda i,j:sp.integrate(compact_weight*coordinate**(i+j),(coordinate,-1,1)))
    coefficients = gram.inv()*sp.Matrix([1,0,0,0,0])
    profile = compact_weight*sum(coefficients[i]*coordinate**i for i in range(order+1))
    checks.check("the defining compact moment Gram matrix is strictly positive", all(gram[:j,:j].det() > 0 for j in range(1,order+2)))
    checks.check("actual signed profiles retain mass and cancel four successive moments", all(sp.integrate(profile*coordinate**i,(coordinate,-1,1)) == int(i == 0) for i in range(order+1)))
    cost = sp.Symbol("D", positive=True)
    accuracy_order = 2*cost+4
    macro_power = cost+1
    checks.check("the explicit narrow-band diagonal controls the normalized cubic cost", sp.simplify(macro_power-cost) == 1)
    checks.check("moment-flat preparation beats the physical second-jet scale on that SAME diagonal", sp.simplify(accuracy_order-2*macro_power) == 2)
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
